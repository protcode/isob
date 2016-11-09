"""
This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2016 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive / Fusion
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This file mainly deals with the handling of the command line data and reading
configuration files converting parameter values as directed.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob

"""

import os
import sys

# python imports
import copy
import datetime
import numpy as np
from pathlib import Path
import socket

# CommonUtils imports
sys.path.insert(0, '..')
from CommonUtils.tools import Stopwatch
from CommonUtils.averagine import AveragineModel
from CommonUtils.hdf5Base import hdf5Base
from CommonUtils.hdf5Mascot import HDF5Mascot
from CommonUtils.hdf5Config import HDF5Config
from CommonUtils.XICprocessor import XIC
from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.IsotopicPatternHandler import IsotopicPatternHandler
from CommonUtils.LoggingManager import Logger
from CommonUtils.MascotModificationsHandler import MascotModifications
import CommonUtils.ExceptionHandler as ExHa
# import CommonUtils.progressbar as progBar
import CommonUtils.QuantMethodHandler as qmh
from CommonUtils.progressReport import progressReport

# project imports
from PCMHandler import PCMdata
from IsotopeCollection import IsotopeCollection
from OutputHandler import OutputHandler

# reportIsos = [-1,0,1,2,3,4,5]


class MS1QuantController:
    def __init__(self, hdf5obj, hdfMascot, cfg=None):
        """
        @brief creates the MS1Q class and links to an open hdf5Path
        @param hdf5Path <mascotHDF5 object>:
        @return:
        """

        self.hdf = hdf5obj
        self.hdfMascot = hdfMascot

        if not cfg:
            self.cfg = ConfigManager('./MS1QuantConfig.cfg')
        else:
            self.cfg = cfg

        self.OutputHandler = OutputHandler(self.hdf, self.cfg, self.hdfMascot.importGroup)

        # self.cfg.scalePpmMda() // has been moved into main()
        # conversionList = ['hyperplexed', 'labelnames']
        #        self.cfg.evalParameters(conversionList)

        self.hdfMascot.checkDatFilePresence(cfg.parameters['runtime']['datfile'])
        self.loadModifications()
        self.modsHandler = MascotModifications(self.allMods)
        #        self.XICprocessor = XICprocessor(hdf5obj, self.cfg)
        self.XICprocessor = XIC(hdf5obj, self.cfg)
        self.IsoPatternHandler = IsotopicPatternHandler(self.hdf, self.cfg, self.XICprocessor)
        self.fsFitting = None
        self.fsSeparate = None
        self.fsTogether = None

        self.xicbins_currentstart = 0
        self.xicbins_nextstart = 0
        self.xicstep = 1000
        self.curr_xic_bins = np.asarray([])
        self.next_xic_bins = np.asarray([])

        self.useSecondary = cfg.parameters['models']['compare']
        if self.cfg.parameters['models']['primary_model'].lower() == 'exact':
            self.primary = 'exact'
            self.secondary = 'averagine'
        else:
            self.primary = 'averagine'
            self.secondary = 'exact'

        self.cfg.parameters['models']['secondary_model'] = self.secondary

    def doMS1quantification(self):
        """
        @brief controls the MS1 quantification methods
        @return:
        """

        self.hdf.createTable(self.hdfMascot.importGroup, 'config_ms1', 'Config', True)
        self.hdfMascot.writeConfig(cfg.convertConfig(), 'config_ms1')
        self.clearQuanTables()
        self.loadData()
        self.processData()

    def clearQuanTables(self):
        ig = self.hdfMascot.importGroup
        self.hdf.removeTable('/%s/quan' % ig)
        self.hdf.removeTable('/%s/quancomp' % ig)
        self.hdf.removeTable('/%s/quanextra' % ig)

    def loadModifications(self):
        """
        @brief loads modification labelDef from Mascot search and pyMSsafe analysis
        @return:
        """
        hdf = self.hdf

        varMods, fixedMods = self.hdfMascot.readModDicts()
        self.fixedMods = fixedMods
        self.allMods = varMods.copy()
        self.allMods.update(fixedMods)

        labellingIsotopeDefs = self.hdfMascot.readIsotopes()
        # labellingIsotopeDefs = hdf.readTable('/rawdata/isotopes')
        self.methID = labellingIsotopeDefs[0]['method_id']
        self.label2id = {}
        for labelDef in labellingIsotopeDefs:
            self.label2id[labelDef['name']] = labelDef['iso_id']
            if 'n-ter' in labelDef['amino'].lower():
                labelDef['amino'] = '+'
            elif 'c-ter' in labelDef['amino'].lower():
                labelDef['amino'] = '-'

        quantMethods = qmh.QuantMethods()
        methodData = quantMethods.getMethodByID(self.methID)

        for labellingMethod in cfg.parameters['general']['labelnames']:
            if self.methID in cfg.parameters['general']['labelnames'][labellingMethod]:
                break

        # filter for modifications
        ms1QuantMods = []
        ms1QuantModIDs = [x for x in varMods if methodData['mascot_name'] in varMods[x]['name']]
        ms1QuantModIDsFixed = [x for x in fixedMods if methodData['mascot_name'] in fixedMods[x]['name']]
        ms1QuantModIDs += ms1QuantModIDsFixed

        if len(ms1QuantModIDsFixed) > 0:
            self.hasFixedLabels = True
        else:
            self.hasFixedLabels = False

        for isotope in labellingIsotopeDefs:
            if isotope['intmz'] == 0:
                # this label has no modification so must be added to the ms1QuantMods
                isoName = '%s:%s+0' % (methodData['mascot_name'], isotope['amino'])
                qm = dict(name=isoName, iso_id=isotope['iso_id'], elem={'C': 0},
                          label=isotope['name'], delta=0, id='0', amino=isotope['amino'])
                ms1QuantMods.append(qm)

        # if methodData['mascot_name'] in ['SILAC8']:
        #     ms1QuantMods.append({'name': 'SILAC: 0', 'iso_id': 35,
        #                          'elem': {'C': 0},
        #                          'label': 'LIGHT', 'delta': 0, 'id': '0', 'amino': 'KR'})
        # elif methodData['mascot_name'] in ['S2ME']:
        #     ms1QuantMods.append({'name': 'S2ME:R+0', 'iso_id': 77,
        #                          'elem': {'C': 0},
        #                          'label': 'S0 DM0', 'delta': 0, 'id': '0', 'amino': 'R'})
        #     ms1QuantMods.append({'name': 'S2ME:R+0', 'iso_id': 78,
        #                          'elem': {'C': 0},
        #                          'label': 'S0 DM6', 'delta': 0, 'id': '0', 'amino': 'R'})

        # combine the labelDef from Mascot with pyMSsafe
        for id in ms1QuantModIDs:
            if id in varMods:
                quantMod = varMods.pop(id)
            else:
                quantMod = fixedMods.pop(id)

            quantMod['id'] = id

            if quantMod['amino'] in ('N-term', 'N_term'):
                quantMod['amino'] = '+'
            elif quantMod['amino'] in ('C-term', 'C_term'):
                quantMod['amino'] = '-'

            for labelDef in labellingIsotopeDefs:
                if labelDef['intmz'] == int(quantMod['delta'] + 0.5) and labelDef['amino'] in quantMod['amino']:
                    quantMod['label'] = labelDef['name']
                    quantMod['iso_id'] = labelDef['iso_id']
                    ms1QuantMods.append(copy.deepcopy(quantMod))
                    # if len(quantMod['amino']) > 1 and not quantMod['amino'] in ('N-term', 'C-term'):
                    if len(quantMod['amino']) > 1:
                        ms1QuantMods[-1]['amino'] = labelDef['amino']

        # keep a copy of the mods list and create a dictionary of the labelDef indexed by the mascot mod_id
        self.ms1QuantMods = ms1QuantMods
        ms1QuantModsIndexed = {}
        for mod in ms1QuantMods:
            if mod['id'] in ms1QuantModsIndexed:
                if not mod['amino'] in ms1QuantModsIndexed[mod['id']]['amino']:
                    ms1QuantModsIndexed[mod['id']]['amino'] += mod['amino']
            else:
                ms1QuantModsIndexed[mod['id']] = copy.deepcopy(mod)

        self.ms1QuantModsIndexed = ms1QuantModsIndexed

        self.varMods = varMods
        self.hyperplexed = methodData['hyperplexed']

        sites = []
        labelDeltas = {}
        for mod in self.ms1QuantMods:
            for aa in mod['amino']:
                if aa not in sites:
                    sites.append(aa)
                labelDeltas.setdefault(mod['label'], {})[aa] = int(mod['delta'] + 0.5)

        # add missing labelDef from pyMssafe - usually just the unlabelled version
        for labelDef in labellingIsotopeDefs:
            if labelDef['name'] in labelDeltas and labelDef['amino'] in labelDeltas[labelDef['name']]:
                continue
            else:
                labelDeltas.setdefault(labelDef['name'], {})[labelDef['amino']] = int(labelDef['mz'] + 0.5)

        self.sites = sites
        self.labelDeltas = labelDeltas

        return

    def loadData(self):
        """
        @brief loads the peptide and MS data
        @return:
        """
        IsoCollect = IsotopeCollection(self.cfg, self.hyperplexed, fixedMods=self.fixedMods,
                                       ms1QuantModifications=self.ms1QuantMods)

        tmpPeptides = self.hdf.getDataEqual('/%s/peptides' % self.hdfMascot.importGroup, 'pepno', 1)
        data = tmpPeptides[:]

        idx = 0
        for tPep in tmpPeptides:
            if tPep['score'] >= self.cfg.parameters['general']['minscore']:
                data[idx] = tPep
                idx += 1

        peptides = data[:idx]
        #a = [d for d in peptides if d['sequence'] in ['HSMNPFCEIAVEEAVR','AQVEIVTDGEEPAEMIQVLGPK','MSAAQALK','EIMLAAKK','GASGTREDPNLVPSISNK','ISMAAALAMNSSPSVR', 'EDPNLVPSISNKR']]
        #a = [d for d in peptides if d['sequence'] == 'LNLKPR']
        #a = [d for d in peptides if d['sequence'] in ('KPATSSKPR', 'LNLKPR')]
  #      a = [d for d in peptides if d['sequence'] in  ['YMEDIQSNVK']]# ,'IIQAHGGTVDPTFTSR','EIAQDFKTDLR']]
    #    a = [d for d in peptides if d['sequence'] in ['IMGLDLPDGGHLTHGYMSDVKR']]#,#'['DDIESLHDKFK','EIAQDFKTDLR','IIQAHGGTVDPTFTSR','YGLNMSSLR','KQSLFAR','VGPLVAGR','KVPELMEMFLPATK','IMGLDLPDGGHLTHGYMSDVKR','AADLLEMATFLQK']]
        #a =  peptides[:200]

        del tmpPeptides
        del data
        #peptides = a


        MSMSheaders = self.hdf.readTable('/rawdata/msmsheader')
        queries = self.hdf.readTable('/%s/queries' % self.hdfMascot.importGroup)

        avFWHM = np.average(MSMSheaders['fwhm'])
        sdFWHM = np.std(MSMSheaders['fwhm'])
        print 'FWHM  = %.3f +/- %.3f' % (avFWHM, sdFWHM)

        # create header location index
        self.headersIndex = {}
        for h in range(len(MSMSheaders)):
            self.headersIndex[MSMSheaders[h]['spec_id']] = h

        # create query to header position index
        self.query2head = {}
        for qry in queries:
            self.query2head[qry['query']] = self.headersIndex[qry['spec_id']]

        # put peptides in PCMdata objects
        allSeqs = set()
        qcSeqs = set()
        allPCMs = {}
        failures = 0
        noMS1quant = 0
        pcmLenghts = {}
        numLabels = {}

        for l in self.labelDeltas:
            numLabels[l] = 0

        progRep = progressReport(len(peptides), self.hdf.filePath.stem, 'loading peptides', 'peptides')
        try:
            # loop through the peptides
            for peptide in peptides:
                progRep.next()
                allSeqs.add(peptide['sequence'])
                pcmlen = len(peptide['sequence'])
                if pcmlen in pcmLenghts:
                    pcmLenghts[pcmlen] += 1
                else:
                    pcmLenghts[pcmlen] = 1

                pcm = PCMdata(peptide, MSMSheaders[self.query2head[peptide['query']]], self.hasFixedLabels)
                self.modsHandler.resolveDualModifiedAminoAcids(pcm)
                pcm.parseMS1QuantMods(self.allMods, self.ms1QuantMods)
                if self.hyperplexed:
                    # this needs hyperplex qc
                    needed = self.hyperplexed
                    qcOK = pcm.doHYPEPLEXqc(self.ms1QuantModsIndexed, self.sites, self.labelDeltas, needed)
                else:
                    # this needs simple QC
                    qcOK = pcm.doMS1QuantQC(self.ms1QuantModsIndexed, self.sites, self.labelDeltas)

                if qcOK:
                    numLabels[pcm.ms1QuantState] += 1
                    if pcm.PCMsf in allPCMs:
                        # the PCM has already been seen
                        allPCMs[pcm.PCMsf].append(pcm)
                    else:
                        allPCMs[pcm.PCMsf] = [pcm]
                        qcSeqs.add(pcm.sequence)

                    IsoCollect.addLookupEntry(pcm)

                elif pcm.ms1QuantState is None:
                    noMS1quant += 1
                else:
                    failures += 1

        except Exception as genEx:
            ExHa.addContext(genEx, 'Error during MS1Quant QC, dealing with peptide %s (Peptide No.: %s)' % (
                peptide.sequence, str(peptide.pepno)))
            raise

        progRep.endReport()

        self.numPeptides = len(peptides)
        self.allSeqs = allSeqs
        self.qcSeqs = qcSeqs
        self.allPCMs = allPCMs
        self.failures = failures
        self.noMS1quant = noMS1quant
        self.pcmLenghts = pcmLenghts
        self.numLabels = numLabels
        self.IsoCollect = IsoCollect

    def processData(self):
        """
        @brief processes the PCM data
        @return:
        """
        msg = 'initialising processing'
        getClustersFromSpectra = self.IsoPatternHandler.getMeasuredClustersFromSpectra
        correctForOverlap = self.IsoPatternHandler.correctForOverlappingIsoPatterns
        getBestFittingIsos = self.IsoPatternHandler.getBestFittingIsoPatternCandidate
        getClustersFromXIC = self.IsoPatternHandler.getMeasuredClustersFromXIC
        useSecondary = self.cfg.parameters['models']['compare']

        # make class variables local
        allPCMs = self.allPCMs
        IsoCollect = self.IsoCollect
        maxXICls = self.cfg.parameters['deisotoping']['maxxicls']
        blindXICacceptLS = self.cfg.parameters['deisotoping']['blindswitch']
        maxSpecLS = self.cfg.parameters['deisotoping']['maxspecls']
        label2id = self.label2id

        if self.methID in self.cfg.parameters['general']['xiconly']:
            calcXIConly = True
        else:
            calcXIConly = False

        print 'Finding XIC cluster data for %i PCMs' % len(allPCMs)
        progRep = progressReport(len(allPCMs), self.hdf.filePath.stem, 'processing XIC', 'precursors')

        # generate ordering by RT for the PCM data
        ordered = []
        #print allPCMs
       # for pcm in allPCMs:
       #     allPCMs[pcm].sort(key=lambda x: x.score, reverse=True)
       #     PCMdata = allPCMs[pcm][0]
       #     ordered.append((PCMdata.RTmsms, pcm))

       # ordered.sort()
        possible_labels = set([x['label'] for x in IsoCollect.ms1QuantModifications])
        self.xicdata = self.hdf.getTable('/rawdata/xicbins')
        self.xicbins_nextstart = 1000

        self.curr_xic_bins = self.xicdata.readWhere('(specid>%s) & (specid<=%s)' % (0, self.xicstep))
        self.next_xic_bins = self.xicdata.readWhere('(specid>%s) & (specid<=%s)' % (self.xicbins_nextstart, self.xicbins_nextstart+self.xicstep))

        #for o in ordered:
        # order all data by best scoring PCMs' Retention time
        #lsfFH=open('LS-fit_data.txt', 'w')
        for PCMdata in sorted([sorted(x, key=lambda y:y.score, reverse=True)[0] for x in allPCMs.values()], key=lambda x: x.RTmsms):
            try:
                progRep.next()
                #PCMdata = allPCMs[o[1]][0]
                #print PCMdata.sequence, PCMdata.PCM, PCMdata.RTmsms
                # generate isotope patterns using exact and averagine models
                cfg.mypcm =  PCMdata.PCM
                msg = 'Generating isotope patterns'
                theoIsoPatterns = IsoCollect.getMS1quantModifiedIsoPatterns(PCMdata)
                bestMethod = None
                prevLabel = None
                associatedIsoPatterns = {'XIC': {},
                                         'priorSurvey': {},
                                         'apexSurvey': {}}
            #    print 'PCMdata.survey', PCMdata.survey, PCMdata.scanApex, PCMdata.survey - PCMdata.scanApex
                if calcXIConly:
                    apexSpec, prevSurvey = None, None
                else:
                    apexSpec, prevSurvey = self.getSpectrumValues(PCMdata.scanApex, PCMdata.survey)

                # read spectrum data for isotope detection
                xicSets = {}
                missingIons = {}
                AS_numpresent = 0
                PS_numpresent = 0

                #oprevSurvey = self.hdf.getDataEqual('/rawdata/xicbins', 'specid', PCMdata.survey)
                #oapexSpec = self.hdf.getDataEqual('/rawdata/xicbins', 'specid', PCMdata.scanApex)

                specIDs = dict(apex=PCMdata.scanApex, prior=PCMdata.survey)

                msg = 'Getting XIC cluster data'
                for label in possible_labels:
                    missingIons[label] = dict(apex=[], prior=[])
                    associatedIsoPatterns['priorSurvey'][label] = {}
                    associatedIsoPatterns['apexSurvey'][label] = {}
                    labelIsotopes = theoIsoPatterns[label]

                    IsoPatternCandidates_XIC = getClustersFromXIC(labelIsotopes, PCMdata)
                    xicSets[label] = copy.deepcopy(IsoPatternCandidates_XIC)
                isoPattern_XIC = self.getOptimalXICSet(xicSets, theoIsoPatterns, PCMdata, maxXICls)

                foundlabels = possible_labels & set([x for x, y in isoPattern_XIC.iteritems() if y])
                XIC_numpresent = len(foundlabels)
                try:
                    if len(possible_labels) == XIC_numpresent:
                        # if we find all possible labels in XIC  found and the fit is good enough to be accepted
                        #  we can stop here
                        if isoPattern_XIC['clusterFitMax'] <= blindXICacceptLS:
                         #   print 'XIC would be it', isoPattern_XIC['clusterFitMax'],
                            bestMethod = 'XIC'

                except KeyError:
                    print 'no lssquares. ', isoPattern_XIC

                secondary_fit_max = 0
                secondary_fit_sum = 0
                for v in isoPattern_XIC:
                    if v in possible_labels:
                        if isoPattern_XIC[v]:
                            if isoPattern_XIC[v]['secondary_fit'] > secondary_fit_max:
                                secondary_fit_max = isoPattern_XIC[v]['secondary_fit']
                            if isoPattern_XIC[v]['secondary_fit'] is not None:
                                secondary_fit_sum += isoPattern_XIC[v]['secondary_fit']
                associatedIsoPatterns['XIC'] = isoPattern_XIC
                associatedIsoPatterns['XIC']['secondaryFitQuality'] = secondary_fit_sum
                associatedIsoPatterns['XIC']['secondaryFitMax'] = secondary_fit_max

                if bestMethod == 'XIC':
                    switched = 0
#                    d = associatedIsoPatterns[bestMethod]
#                    for label in ('HEAVY', 'LIGHT'):
#                        try:
 #                           lsfFH.write('%s\t%s\t%s' % (self.hdf.filePath.stem, PCMdata.spec_id, ''.join(x for x in d[label]['fitdata'])))
  #                      except KeyError:
  #                          pass
                            #print label, 'no fit data'
                    self.IsoPatternHandler.getMisincorporationRatios(associatedIsoPatterns, theoIsoPatterns)
                    self.OutputHandler.PCMoutput(PCMdata, associatedIsoPatterns, bestMethod, theoIsoPatterns, switched,
                                                 missingIons, specIDs, label2id, self.fsFitting, self.fsSeparate,
                                                 self.fsTogether)
                    self.OutputHandler.writeIsotopeData(PCMdata.spec_id, theoIsoPatterns, label2id, associatedIsoPatterns)
                    continue
                else:

                    missingIons = {}
                    sum_AS_LS_secondary = 0
                    max_AS_LS_secondary = 0
                    sum_PS_LS_secondary = 0
                    max_PS_LS_secondary = 0

                    max_PS_LS = 0
                    sum_PS_LS = 0
                    max_AS_LS = 0
                    sum_AS_LS = 0

                    PS_all_LS = [0]
                    AS_all_LS = [0]
                    AS_all_LS_secondary = [0]
                    PS_all_LS_secondary = [0]

                    # extract XIC data and identify isotope patterns also for spectrum data
                    prevLabel = None
                    for label in theoIsoPatterns['byIncreasingMZ']:
                        msg = 'processing label:', label
                        if label in theoIsoPatterns:
                            labelIsotopes = theoIsoPatterns[label]
#                            print 'labelIsotopes', labelIsotopes
                            try:
                                # generate XIC data clusters
                                if calcXIConly:
                                    associatedIsoPatterns['priorSurvey'][label] = {}
                                    associatedIsoPatterns['apexSurvey'][label] = {}
                                    missingIons[label] = dict(apex=[], prior=[])
                                else:
                                    # generate Isotope clusters from spectra
                                    msg = 'getting specta isotopes'
                                    IsoPatternList_PS, PSmissing = getClustersFromSpectra(prevSurvey, labelIsotopes)
                                    IsoPatternList_AS, ASmissing = getClustersFromSpectra(apexSpec, labelIsotopes)
                                    missingIons[label] = dict(apex=ASmissing, prior=PSmissing)
                                    # do isotope correction for spectrum data
                                    msg = 'doing isotope correction for spectrum'
                                    correctForOverlap(IsoPatternList_PS, theoIsoPatterns, label,
                                                                          prevLabel, associatedIsoPatterns['priorSurvey'],
                                                                          1.0)
                                    correctForOverlap(IsoPatternList_AS, theoIsoPatterns, label,
                                                                          prevLabel, associatedIsoPatterns['apexSurvey'],
                                                                          1.0)
                                    # extract the best cluster (should be only one)
                                    msg = 'extracting the best cluster for spectrum'

                                    cfg.model = 'PS'
                                    isoPattern_PS = getBestFittingIsos(IsoPatternList_PS, labelIsotopes, PCMdata.charge,
                                                                       useSecondary)
                                    cfg.model = 'AS'
                                    isoPattern_AS = getBestFittingIsos(IsoPatternList_AS, labelIsotopes, PCMdata.charge,
                                                                       useSecondary)

                                    # append the data to associatedIsoPatterns if the fit is good enough
                                    if isoPattern_PS and (maxSpecLS == 0 or isoPattern_PS['leastSquares'] <= maxSpecLS):
                                        associatedIsoPatterns['priorSurvey'][label] = isoPattern_PS.copy()
                                        associatedIsoPatterns['priorSurvey'][label]['ms1SpecID'] = PCMdata.survey
                                        PS_numpresent += 1
                                        PS_all_LS.append(isoPattern_PS['leastSquares'])
                                        sum_PS_LS += isoPattern_PS['leastSquares']
                                        if isoPattern_PS['leastSquares'] > max_PS_LS:
                                            max_PS_LS = isoPattern_PS['leastSquares']

                                        PS_all_LS_secondary.append(isoPattern_PS['secondary_fit'])
                                    else:
                                        associatedIsoPatterns['priorSurvey'][label] = {}
                                    #print 'isoPattern_AS', isoPattern_AS,
                                   # if 'leastSquares' in  isoPattern_AS:
                                    #        print isoPattern_AS['leastSquares']
                                    if isoPattern_AS and (maxSpecLS == 0 or isoPattern_AS['leastSquares'] <= maxSpecLS):
                                        associatedIsoPatterns['apexSurvey'][label] = isoPattern_AS.copy()
                                        associatedIsoPatterns['apexSurvey'][label]['ms1SpecID'] = PCMdata.scanApex
                                        AS_numpresent += 1
                                     #   print 'AS numprsent plus 1'
                                        AS_all_LS.append(isoPattern_AS['leastSquares'])
                                        sum_AS_LS += isoPattern_AS['leastSquares']
                                        if isoPattern_AS['leastSquares'] > max_AS_LS:
                                            max_AS_LS = isoPattern_AS['leastSquares']
                                        AS_all_LS_secondary.append(isoPattern_AS['secondary_fit'])
                                    else:
                                        associatedIsoPatterns['apexSurvey'][label] = {}

                                prevLabel = label

                            except Exception as genEx:
                                ExHa.addContext(genEx, 'Error during processing of PCM %s, label %s' % (PCMdata.PCM, label))
                                raise

                    # process the XIC data (needs to run on all clusters in all labels to link clusters properly
                    # check whether the misincorporation peaks (if present) of the selected clusters have a reasonable
                    # intensity or indicate an interference with another pattern
                    self.IsoPatternHandler.getMisincorporationRatios(associatedIsoPatterns, theoIsoPatterns)

                    switched = 0

                    maxNumPresent = max(XIC_numpresent, AS_numpresent, PS_numpresent)
                    if self.useSecondary:
                        associatedIsoPatterns['apexSurvey']['secondaryFitQuality'] = sum(AS_all_LS_secondary)
                        associatedIsoPatterns['priorSurvey']['secondaryFitQuality'] = sum(PS_all_LS_secondary)
                        associatedIsoPatterns['apexSurvey']['secondaryFitMax'] = max(AS_all_LS_secondary)
                        associatedIsoPatterns['priorSurvey']['secondaryFitMax'] = max(PS_all_LS_secondary)

                    associatedIsoPatterns['apexSurvey']['patternSetFitQuality'] = sum_AS_LS
                    associatedIsoPatterns['apexSurvey']['clusterFitMax'] = max_AS_LS
                    associatedIsoPatterns['priorSurvey']['patternSetFitQuality'] = sum_PS_LS
                    associatedIsoPatterns['priorSurvey']['clusterFitMax'] = max_PS_LS
                    try:
                        refVal = associatedIsoPatterns['XIC']['patternSetFitQuality']
                    except KeyError:

                        refVal = 0

                    if calcXIConly:
                        bestMethod = 'XIC'
                    elif maxNumPresent > 0:
                        if XIC_numpresent == maxNumPresent:
                            # we could use the XIC
                            if associatedIsoPatterns['XIC']['clusterFitMax'] <= blindXICacceptLS:
                                # this should be used as good fit
                                bestMethod = 'XIC'
                            #    print 'we get to this point if label number is less than optimal but fit still good',associatedIsoPatterns['XIC']['clusterFitMax'], maxNumPresent, foundlabels
                            elif (AS_numpresent == maxNumPresent and
                                  associatedIsoPatterns['apexSurvey']['patternSetFitQuality'] < refVal):
                                # AS has same number of theoIsoPatterns and better fitting
                                bestMethod = 'apexSurvey'
                                switched = 1
                            #    print 'choosing Apex as ASnumpresente and maxnum prsent identical', refVal, AS_numpresent, maxNumPresent
                            elif (PS_numpresent == maxNumPresent and
                                  associatedIsoPatterns['priorSurvey']['patternSetFitQuality'] < refVal):
                                bestMethod = 'priorSurvey'
                             #   print 'choosing ps as best method',AS_numpresent, maxNumPresent
                                switched = 1
                            else:
                                bestMethod = 'XIC'
                              #  print 'choosing XIC as best method; xic count same for all and no other mathces',AS_numpresent, maxNumPresent, PS_numpresent
                        elif AS_numpresent == maxNumPresent:
                            bestMethod = 'apexSurvey'
                            #print 'apexSurvey cos numbers the same and but not XIC count',AS_numpresent, maxNumPresent
                        else:
                            bestMethod = 'priorSurvey'
                            #print 'priorSurvey cos numbers the same and but no XIC count, no apex'
                    else:
                        bestMethod = 'XIC'

                    msg = 'outputting data'

                    # d = associatedIsoPatterns[bestMethod]
                    # for label in ('HEAVY', 'LIGHT'):
                    #
                    #     try:
                    #         lsfFH.write('%s\t%s\t%s' % (self.hdf.filePath.stem,PCMdata.spec_id, ''.join(x for x in d[label]['fitdata'])))
                    #     except KeyError:
                    #         pass
                    #         #print label, 'no fit data'
                    #
                    self.OutputHandler.PCMoutput(PCMdata, associatedIsoPatterns, bestMethod, theoIsoPatterns, switched,
                                                 missingIons, specIDs, label2id, self.fsFitting, self.fsSeparate,
                                                 self.fsTogether)
                    self.OutputHandler.writeIsotopeData(PCMdata.spec_id, theoIsoPatterns, label2id, associatedIsoPatterns)

            except Exception as genEx:

                ExHa.addContext(genEx, 'Error during processing of PCM %s, msg = %s' % (PCMdata.PCM, msg))
                raise

        progRep.endReport()
        #lsfFH.close()
        return

    def getOptimalXICSet(self, xicSets, theoIsoPattern, PCMdata, maxLeastSq):
        """
        @brief Select the best set (i.e., representing the different labels) of isotope clusters detected from XICs.
                The 12C peaks of the clusters are required to overlap. In case of ambiguity, choose the set with most
                labels and lowest summed least squares values.
        @param xicSets: <dictionary> containing the all possible clusters per label
        @param theoIsoPattern: <dictionary> containing the corresponding theoretical isotope patterns
        @param PCMdata: <int> charge state of the corresponding pcm
        @param maxLeastSq: <float>: threshold for maximum accepted least squares value (currently turned off by default)
        @return:
        """

        correctForOverlap = self.IsoPatternHandler.correctForOverlappingIsoPatterns
        getBestFittingIsos = self.IsoPatternHandler.getBestFittingIsoPatternCandidate
        minSwitch = cfg.parameters['deisotoping']['blindswitch']

        presentLabels = []
        for lbl in theoIsoPattern['byIncreasingMZ']:
            if xicSets[lbl]:
                presentLabels.append(lbl)

        SetCollection = []
        for pL in presentLabels:
            extras = []
            for cluster in xicSets[pL]:
                #print 'cluster', cluster
                if pL == presentLabels[0]:
                    # first label: so all clusters turned into starting sets
                    newSet = {pL: cluster, 'refRTleft': cluster[0]['rt50left'],
                              'refRTright': cluster[0]['rt50right'], 'numLabels': 1}
                    SetCollection.append(newSet)

                else:
                    missing = 1
                    for set in SetCollection:
                        if self.compatibleRTs(set, cluster[0]):
                            # found compatible Cluster
                            missing = 0
                            if pL in set:
                                # has existing match so duplicate set
                                extras.append(copy.deepcopy(set))
                                extras[-1][pL] = copy.deepcopy(cluster)
                            else:
                                # append new label data to set
                                set[pL] = copy.deepcopy(cluster)
                                set['numLabels'] += 1

                    if missing:
                        # orphaned clusters should start a new group
                        newSet = {pL: cluster, 'refRTleft': cluster[0]['rt50left'],
                                  'refRTright': cluster[0]['rt50right'], 'numLabels': 1}
                        SetCollection.append(newSet)

            if extras:
                # make sure any additional sets are incorporated in the main list
                SetCollection.extend(extras)

        if not SetCollection:
            # print 'empty SetCollection'
            # print 'XIC sets', str(xicSets)
            empty = {}
            for iso in theoIsoPattern['byIncreasingMZ']:
                empty[iso] = {}
            return empty

        # SetCollection.sort(key=lambda x: x['numLabels'], reverse=True)

        # maxMatched = SetCollection[0]['numLabels']
        #        XICclusters = []
        for set in SetCollection:
            #            if set['numLabels'] == maxMatched:
            del set['numLabels']
            del set['refRTleft']
            del set['refRTright']
            #            XICclusters.append(set)

        # need to apply overlap correction and calculate least squares fit
        for set in SetCollection:
            prevLabel = None
            for label in presentLabels:
                if prevLabel and prevLabel in set and label in set and set[prevLabel]:
                    overlapInfluence = self.calculateOverlapInfluence(set[label][0], set[prevLabel]['cluster'][0])
                    # overlapInfluence = 1.0
                else:
                    overlapInfluence = 1.0

                if label in set:
                    set[label]['overlapinfluence'] = overlapInfluence
                    labelCopy = copy.deepcopy(set[label])
                    overlapped = correctForOverlap([set[label]], theoIsoPattern, label, prevLabel, set,
                                                   overlapInfluence)

                    cfg.thislab = label
                    cfg.model = 'XIC'
                    bestXICcluster = getBestFittingIsos(overlapped, theoIsoPattern[label], PCMdata.charge,
                                                        self.useSecondary)
                    # comparison to full overlap data
                    if bestXICcluster and overlapInfluence < 1 and bestXICcluster['leastSquares'] > minSwitch:
                        over2 = correctForOverlap([labelCopy], theoIsoPattern, label, prevLabel, set, 1.0)
                        cfg.model = 'XIC2'
                        bestXIC2 = getBestFittingIsos(over2, theoIsoPattern[label], PCMdata.charge, self.useSecondary)

                        if bestXIC2 and bestXIC2['leastSquares'] < bestXICcluster['leastSquares']:
                            bestXICcluster = copy.deepcopy(bestXIC2)

                    if bestXICcluster and (bestXICcluster['leastSquares'] <= maxLeastSq or maxLeastSq == 0):
                        set[label] = bestXICcluster
                    else:
                        set[label] = {}

                prevLabel = label

        for set in SetCollection:
            nonEmptys = 0
            for k in theoIsoPattern['byIncreasingMZ']:
                if set.setdefault(k, {}):
                    nonEmptys += 1
            set['nonEmptys'] = nonEmptys

        SetCollection.sort(key=lambda x: x['nonEmptys'], reverse=True)
        maxNonEmptys = SetCollection[0]['nonEmptys']

        best = SetCollection.pop(0)
        self.getPatternSetFitQuality(best, presentLabels)

        #        numChoices = 1
        for set in SetCollection:
            if set['nonEmptys'] == maxNonEmptys:
                self.getPatternSetFitQuality(set, presentLabels)
                if set['patternSetFitQuality'] < best['patternSetFitQuality']:
                    best = set

        return best

    def calculateOverlapInfluence(self, current12C, prev12C):

        if current12C['rt'] > prev12C['rt']:
            base = 2 * (prev12C['rt50right'] - prev12C['rt'])
            deltaApex = current12C['rt'] - prev12C['rt']
        else:
            base = 2 * (prev12C['rt'] - prev12C['rt50left'])
            deltaApex = prev12C['rt'] - current12C['rt']

        if deltaApex > base:
            ovelapInfluence = 0
        else:
            ovelapInfluence = (base - deltaApex) / base
        return ovelapInfluence

    def getPatternSetFitQuality(self, cluster, labels):
        """
        @brief calculate the sum and maximum leastSquares values for the clusters
        @param cluster <dictionary>: containing the set of cluster data for all lables
        @param labels <list>: containing the labels present in the cluster data
        @return:
        """
        summedLS = 0.0
        maxLS = 0
        for label in labels:
            if label in cluster and cluster[label]:
                summedLS += cluster[label]['leastSquares']
                if cluster[label]['leastSquares'] > maxLS:
                    maxLS = cluster[label]['leastSquares']

        cluster['patternSetFitQuality'] = summedLS
        cluster['clusterFitMax'] = maxLS

    def compatibleRTs(self, ref, test):
        allowLeft = self.cfg.parameters['xic']['rtallowleft']
        allowRight = self.cfg.parameters['xic']['rtallowright']

        if (ref['refRTleft'] - allowLeft) <= test['rt50right'] and (ref['refRTright'] + allowRight) >= test['rt50left']:
            return True
        else:
            return False

    def getSpectrumValues(self, apexspecid, surveyspecid):
        movetonext = 0
        if apexspecid < self.xicbins_currentstart + self.xicstep:
            q = self.curr_xic_bins[self.curr_xic_bins['specid'] == apexspecid]

            #print 'in apex current', self.xicbins_currentstart,'-->',self.xicbins_currentstart + self.xicstep
        elif apexspecid < self.xicbins_nextstart + self.xicstep:
            q = self.next_xic_bins[self.next_xic_bins['specid'] == apexspecid]
            movetonext += 1
           # print 'in apex next',self.xicbins_nextstart,'-->', self.xicbins_nextstart + self.xicstep
        else:
          #  print 'apex not found'
            q = self.xicdata.readWhere('(specid==%s)' % apexspecid)
            movetonext += 1
        #print 'apexspecid', apexspecid, self.xicbins_nextstart + self.xicstep, self.xicbins_currentstart + self.xicstep
        if surveyspecid < self.xicbins_currentstart + self.xicstep:
            p = self.curr_xic_bins[self.curr_xic_bins['specid'] == surveyspecid]
         #   print 'in survey current', self.xicbins_currentstart, '-->', self.xicbins_currentstart + self.xicstep
        elif surveyspecid < self.xicbins_nextstart + self.xicstep:
            p = self.next_xic_bins[self.next_xic_bins['specid'] == surveyspecid]
            #print 'in survey next', self.xicbins_nextstart,'-->', self.xicbins_nextstart + self.xicstep
            movetonext += 1
        else:
            p = self.xicdata.readWhere('(specid==%s)' % surveyspecid)
            movetonext+=1
           # print 'in survey not found', self.xicbins_nextstart + self.xicstep
        #print 'movetonext ', movetonext
        if movetonext == 2:
            #print 'change'
            sys.stdout.flush()
            self.curr_xic_bins = self.next_xic_bins
            self.xicbins_currentstart = self.xicbins_nextstart
            self.xicbins_nextstart += self.xicstep
            self.next_xic_bins = self.xicdata.readWhere('(specid>%s) & (specid<=%s)' %
                                                         (self.xicbins_nextstart, self.xicbins_nextstart+self.xicstep))
        if not len(q):
            q = self.xicdata.readWhere('(specid==%s)' % apexspecid)
            #print 'eh why didn''t I find q data first time (apex) %s' % apexspecid
        if not len(p):
            p = self.xicdata.readWhere('(specid==%s)' % surveyspecid)
           # print 'eh why didn''t I find p data first time (surveyspecid) %s' % surveyspecid
        return q, p



if __name__ == '__main__':

    logs = 0
    appVers = 0
    hdf5Path = ''
    try:
        cfg = ConfigManager('./MS1QuantConfig.cfg')
        configErr = cfg.evaluateCommandLineArgs(sys.argv)

        cfg.scalePpmMda()

        logParam = cfg.parameters['logging']
        logPath = logParam['logdir'].joinpath(cfg.parameters['general']['pid'])
        if not logPath.exists():
            logPath.mkdir(parents=True)
        logFile = logPath.joinpath(logParam['logfile'])
        logs = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], False)
        logs.setLog('MS1quant')

        hdf5Path = cfg.parameters['runtime']['datadir'].joinpath(cfg.parameters['runtime']['filename'])

        hdf5File = hdf5Base(hdf5Path, [])
        hdf5File.appendOpen()

        hdfMascot = HDF5Mascot(hdfBaseObj=hdf5File)
        importGroup = hdfMascot.checkDatFilePresence(cfg.parameters['runtime']['datfile'])

        if not importGroup:
            raise ExHa.HDFmissingDataError('DAT file (%s) missing from HDF5 file (%s)' %
                                           (cfg.parameters['runtime']['datfile'],
                                            cfg.parameters['runtime']['filename']))

        hdfMascot.getMSconfig(cfg)

        # if hdf5 file already contains quant data, remove it
        for quanTable in ['quan', 'quanextra', 'quaniso', 'quanerror']:
            # hdf5File.removeTable('/rawdata/' + quanTable)
            hdf5File.removeTable('/%s/%s' % (importGroup, quanTable))

        sw_total = Stopwatch()

        MS1Q = MS1QuantController(hdf5File, hdfMascot, cfg)

        # output key parameters
        print 'Processing:           %s' % cfg.parameters['runtime']['filename']
        print 'From:                 %s' % cfg.parameters['runtime']['datadir']
        print 'Minimum Mascot score: %.1f' % MS1Q.cfg.parameters['general']['minscore']
        print 'Switch threshold:     %.2f' % MS1Q.cfg.parameters['deisotoping']['blindswitch']

        if cfg.parameters['general']['suffix'] == 'not_set':
            cfg.parameters['general']['suffix'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # if output:
        if cfg.parameters['general']['txtoutput']:
            print 'Textfile output is ON'
            print 'extension:            "%s"' % cfg.parameters['general']['suffix']
            MS1Q.OutputHandler.createOutputFiles(cfg.parameters['general']['suffix'])

        MS1Q.doMS1quantification()
        hdf5File.close()

        MS1Q.OutputHandler.showSummary(MS1Q.allPCMs, MS1Q.pcmLenghts, MS1Q.numLabels, MS1Q.IsoCollect, MS1Q.numPeptides,
                                       MS1Q.allSeqs, MS1Q.qcSeqs, MS1Q.failures, MS1Q.noMS1quant)

        print 'finished'
        sw_total.stop()
        sw_total.write()

        if cfg.parameters['general']['txtoutput']:
            MS1Q.OutputHandler.closeFiles()

    except Exception as genEx:
        ExHa.reformatException(genEx)
        ExHa.addContext(genEx, 'Error during MS1 quantification')
        ExHa.exportError2File(genEx, cfg.parameters['runtime']['datadir'] / 'errors.txt')
        if logs:
            logs.log.warning(ExHa.oneLineRepr(genEx))
        else:
            print ExHa.multiLineRepr(genEx)
