"""
This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This file mainly facilitates the storage and manipulation of data read from
Mascot .dat files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/cellzome/isobarquant
"""

import sys
from xml.sax import make_parser
import numpy as np
import re

# CommonUtils imports
sys.path.insert(0, '..')

from CommonUtils.MathTools import Statistics
from CommonUtils.tools import *
import CommonUtils.ExceptionHandler as ExHa

# project imports
import peptide
from unimodxmlparser import UnimodXMLParser

rx_pepmaindata = re.compile('q[0-9]+_p[0-9]+$')


class Datfile:
    '''
    Class representing the information from a Mascot datfile
    '''

    def __init__(self, filedic, hdf5, cfg, logs, searchID, quantMethod):
        '''
        @brief initialise the Datfile object
        @param filedic <dictionary>: containig the full path to the dat file and other file data
        @param hdf5 <hdf5 object>: containing all the writing methods for hdf5 output
        @param cfg <cfg object>: containing all the running parameters
        @param dbobj <db object>: linking to the meet database
        '''
        self.filedic = filedic
        self.searchid = searchID
        self.quantMethod = quantMethod
        datfilename = filedic['dat']
        if datfilename == str(None):
            raise ExHa.FileNotFoundException('no resulting datfile found with search')
        self.datfilename = datfilename
        self.cfg = cfg
        self.logs = logs
        self.hdfMascot = hdf5
        self.dataDir = filedic['datpath'].parent

        self.spectra = {}
        self.peptidecounter = 0
        self.failedsequences = 0
        self.sequences = {}
        self.seq2acc = {}
        self.hookpeps = []
        self.hookppm = []

        self.analtimes = dict(first=40, last=0, early=[])
        self.stats = dict(numpeps=0, numfailedpeps=0, numspectra_nopeps=0)

    def addMasses(self, massdict):
        '''
        @brief takes massdict and converts data to numbers, also groups modification data together, saves data to
               HDF5 file
        @param massdict <dictionary>: containing parsed data from masses section
        '''
        # todo: put the ptm's into a separate dictionary
        highIndexMods = {}
        for idx in range(10, 33):
            highIndexMods[str(idx)] = chr(idx + 87)

        keys = massdict.keys()
        mods = [x for x in keys if x[:5] == 'delta']
        fixed = [x for x in keys if x[:8] == 'FixedMod' and x[:9] != 'FixedModR']

        # convert all data to floats where possible
        for key in massdict:
            try:
                massdict[key] = float(massdict[key])
            except:
                continue

        moddict = {}

        for m in mods:
            modnum = m[5:]
            if modnum in highIndexMods:
                modnum = highIndexMods[modnum]
            mod = dict(da_delta=0.0, name='', amino='', nl=0.0, nlm=0.0, modtype='variable')
            data = massdict[m].split(',')
            mod['da_delta'] = float(data[0])
            mod['id'] = modnum
            end = data[1].rfind('(')
            if end > 0:
                mod['name'] = data[1][:end - 1]
                end2 = data[1].find(')', end)
                mod['amino'] = data[1][end + 1: end2]
            else:
                # catch for iTRAQ mods held in braces!
                mod['name'] = data[1][1:-1]
                mod['amino'] = data[1][1]

            # set the relevant mod status
            mod['relevant'] = self.isRelevantMod(mod['name'])

            # now extract NeutralLoss data
            nl = 'NeutralLoss' + modnum
            if nl in massdict:
                mod['nl'] = float(massdict[nl])
                del massdict[nl]
            nl += '_master'
            if nl in massdict:
                mod['nlm'] = float(massdict[nl])
                del massdict[nl]

            del massdict[m]
            moddict[modnum] = mod.copy()

        # handle the fixed modification
        for f in fixed:
            modnum = f[8:]
            mod = dict(da_delta=0.0, name='', amino='', nl=0.0, nlm=0.0, modtype='fixed')
            data = massdict[f].split(',')
            mod['da_delta'] = float(data[0])
            end = data[1].find('(')
            mod['name'] = data[1][:end - 1]
            del massdict[f]

            # set the relevant mod status: all fixed mods are non-relevant
            mod['relevant'] = 0

            res = 'FixedModResidues' + modnum
            if res in massdict:
                mod['amino'] = massdict[res]
                del massdict[res]

            if mod['amino'] == 'N_term':
                # N-terminal modification
                mod['id'] = '+'
                moddict['+'] = mod.copy()
            elif mod['amino'] == 'C_term':
                # C-terminal modification
                mod['id'] = '-'
                moddict['-'] = mod.copy()
            else:
                mod['id'] = mod['amino']
                moddict[mod['amino']] = mod.copy()

        massdict['Proton'] = 1.007276

        # store the masses in the hdf5 file
        self.hdfMascot.writeMasses(massdict)
        if moddict:
            self.hdfMascot.writeMods(moddict)

        self.masses = massdict
        self.mods = moddict

        self.logs.datlog.debug('masses are %s ' % str(massdict))

    def isRelevantMod(self, name):
        '''
        @brief tests if the modification is a relevant one or not
        @param name <string>: the name of the testing modification
        '''
        # regex searching expects uppercase so force conversion
        lcName = name.lower()
        if self.quantMethod:
            quan = self.quantMethod['meth_type'].lower()
        else:
            quan = 'xxx'

        if quan in lcName:
            return 0
        elif lcName in self.cfg.parameters['modifications']['ignoremods']:
            return 0
        else:
            return 1

    def addUnimod(self, unimod):
        '''
        @brief saves the unimod xml data into a temp file then parses the data
        @param unimod <list>: containing all the xml data for parsing
        '''

        tmpXML = self.dataDir.joinpath(self.cfg.parameters['general']['tmpxml'])
        xmlOut = open(str(tmpXML), 'w')

        # strip leading & trailing non xml lines
        while unimod[0][0] != '<':
            unimod.pop(0)
        while unimod[-1][0] != '<':
            unimod.pop(-1)
        for x in unimod:
            xmlOut.write(x)
        xmlOut.close()

        xmlIn = file(str(tmpXML), 'r')
        xml = UnimodXMLParser()
        saxparser = make_parser()
        saxparser.setContentHandler(xml)
        saxparser.parse(xmlIn)
        xmlIn.close()
        tmpXML.unlink()

        self.hdfMascot.writeUnimodAminos(xml.aminoacids)
        self.hdfMascot.writeUnimodElements(xml.elements)

        if hasattr(xml, 'mods'):
            self.hdfMascot.writeUnimodModifications(xml.mods)

            # link the unimod data with the mascot modification data
            self.hdfMascot.updateModsWithElements(xml.mods)

    def addSummary(self, summary):
        '''
        @brief organises the summary data
        @params summary <dictionary>: containint the parsed summary data
        '''
        # find all qmass data lines
        specs = [(int(x[5:]), x[5:], x) for x in summary.keys() if x[:5] == 'qmass']
        specs.sort()
        self.stats['numspectra'] = len(specs)

        # update log and db
        self.logs.datlog.info('Parsing %d summary spectra' % (len(specs)))

        # process data
        for query, querystr, sp in specs:
            spdata = dict(query=query, prec_neutmass=float(summary[sp]), numpeps=0)

            # expand to the qexp line for the same spectrum
            if 'qexp' + querystr in summary:
                expdata = summary['qexp' + querystr].split(',')
                spdata['prec_mz'] = float(expdata[0])
                spdata['prec_charge'] = int(expdata[1].replace('+', ''))

            # add number of db matches
            if 'qmatch' + querystr in summary:
                spdata['matches'] = int(summary['qmatch' + querystr])

            # add homology score threshold
            if 'qplughole' + querystr in summary:
                spdata['homology'] = float(summary['qplughole' + querystr])

            # retain info
            self.spectra[query] = spdata

    def addETsummary(self, summary):
        '''
        @brief organises the summary data
        @params summary <dictionary>: containint the parsed summary data
        '''
        # find all qmass data lines
        specs = [(int(x[5:]), x[5:], x) for x in summary.keys() if x[:5] == 'qmass']
        specs.sort()
        self.stats['numspectra'] = len(specs)

        # update log and db
        self.logs.datlog.info('Parsing %d summary spectra' % (len(specs)))

        # process data
        for query, querystr, sp in specs:
            spdata = dict(query=query, prec_neutmass=float(summary[sp]), numpeps=0)

            # expand to the qexp line for the same spectrum
            if 'qexp' + querystr in summary:
                expdata = summary['qexp' + querystr].split(',')
                spdata['prec_mz'] = float(expdata[0])
                spdata['prec_charge'] = int(expdata[1].replace('+', ''))

            # add number of db matches
            if 'qmatch' + querystr in summary:
                spdata['matches'] = int(summary['qmatch' + querystr])

            # add homology score threshold
            if 'qplughole' + querystr in summary:
                spdata['homology'] = float(summary['qplughole' + querystr])

            # retain info
            self.spectra[query] = spdata

    def addCalPeptide(self, peps4spec):
        '''
        @brief processes the peptide data for a single spectrum
        @param peps4spec <dictionary>: containing all peptides for one spectrum
        '''
        stats = self.stats
        cfg = self.cfg
        pepkeys = [x for x in peps4spec.keys() if rx_pepmaindata.search(x)]
        pepkeys.sort()

        bits = pepkeys[0].split('_')
        query = int(bits[0][1:])

        spec = self.spectra[query]

        self.logs.datlog.debug('Query %i has %i peptides' % (query, len(pepkeys)))

        if len(pepkeys) > 1 and pepkeys[1][-2:] == '10':
            key = pepkeys.pop(1)
            pepkeys.append(key)

        for pep in pepkeys:
            # collect all the data for each peptide and create Peptide objects
            rx_pepdata = re.compile('^' + pep + '(?![0-9])')
            allpeptidedata = {}
            for x in peps4spec.keys():
                if rx_pepdata.search(x):
                    allpeptidedata[x] = peps4spec[x]

            # first create the peptide instance
            pepobj = peptide.Peptide(allpeptidedata, self)

            # then test if the sequence is valid
            if pepobj.isValidSequence():
                # only include peptides if the sequence is valid
                self.spectra[query]['numpeps'] += 1
                stats['numpeps'] += 1

                for attr in ['sequence', 'pepno', 'modsVariable', 'modsFixed', 'mass', 'da_delta', 'score']:
                    spec[attr] = getattr(pepobj, attr)

                accessions = [x['accession'] for x in pepobj.proteins]

                if len([x for x in accessions if x.startswith('DD') or x.startswith('###REV###') or
                        x.startswith('###RND###')]) == len(accessions):
                    hitType = 'REV'
                else:
                    hitType = 'FWD'
                spec['hitType'] = hitType

                # only process one peptide per query
                break
            else:
                stats['numfailedpeps'] += 1

    def addPeptides(self, peps4spec):
        '''
        @brief processes the peptide data for a single spectrum
        @param peps4spec <dictionary>: containing all peptides for one spectrum
        '''
        stats = self.stats
        cfg = self.cfg
        pepkeys = [x for x in peps4spec.keys() if rx_pepmaindata.search(x)]
        pepkeys.sort()

        bits = pepkeys[0].split('_')
        query = int(bits[0][1:])
        self.logs.datlog.debug('Query %i has %i peptides' % (query, len(pepkeys)))

        if len(pepkeys) > 1 and pepkeys[1][-2:] == '10':
            key = pepkeys.pop(1)
            pepkeys.append(key)
        peplist = []
        seq2acc = self.seq2acc
        sequences = self.sequences
        self.peptidecounter += len(pepkeys)

        for pep in pepkeys:
            # collect all the data for each peptide and create Peptide objects
            rx_pepdata = re.compile('^' + pep + '(?![0-9])')
            allpeptidedata = {}
            for x in peps4spec.keys():
                if rx_pepdata.search(x):
                    allpeptidedata[x] = peps4spec[x]

            # first create the peptide instance
            pepobj = peptide.Peptide(allpeptidedata, self)

            # then test if the sequence is valid
            if pepobj.isValidSequence():
                # only include peptides if the sequence is valid
                self.spectra[query]['numpeps'] += 1
                stats['numpeps'] += 1
                peplist.append(pepobj)
            else:
                stats['numfailedpeps'] += 1

        # now do the QC of the peptide set
        if peplist:
            self.doPeptideSetQC(peplist)

        hasHook = ''
        for pep in peplist:
            if pep.useinprot == 0 or pep.retain == 0:
                continue
            seq = pep.sequence
            score = pep.score
            if pep.is_hook:
                hookscore = score
                hasHook = ', has hook peptide'
            else:
                hookscore = 0.0

            # build dictionary of pep sequence to protein accession
            if seq in sequences:
                sequences[seq] += 1
                # accumulate data independently
                seq2acc[seq]['numpep'] += 1
                if pep.is_hook > seq2acc[seq]['hook']:
                    seq2acc[seq]['hook'] = pep.is_hook

                if hookscore > seq2acc[seq]['hookscore']:
                    seq2acc[seq]['hookscore'] = hookscore

                if score > seq2acc[seq]['pepscore']:
                    seq2acc[seq]['pepscore'] = score

                if pep.pepno < seq2acc[seq]['bestczrank']:
                    seq2acc[seq]['bestczrank'] = pep.pepno
            else:
                sequences[seq] = 1
                seq2acc[seq] = dict(prots=pep.proteins[:], hook=pep.is_hook, numpep=1,
                                    hookscore=hookscore, pepscore=score, bestczrank=pep.pepno)

        self.logs.datlog.debug('%i peptides pass QC%s' % (len(peplist), hasHook))

    def createETpeptidesTable(self):
        '''
        @brief creates teh ETpeptides table in the HDF5 file
        '''
        self.hdfMascot.createETpeptidesTable()

    def addETpeptides(self, peps4spec):
        '''
        @brief processes the peptide data for a single spectrum
        @param peps4spec <dictionary>: containing all peptides for one spectrum
        '''
        stats = self.stats

        pepkeys = [x for x in peps4spec.keys() if rx_pepmaindata.search(x)]
        pepkeys.sort()

        bits = pepkeys[0].split('_')
        query = int(bits[0][1:])

        if len(pepkeys) > 1 and pepkeys[1][-2:] == '10':
            key = pepkeys.pop(1)
            pepkeys.append(key)
        peplist = []
        seq2acc = self.seq2acc
        sequences = self.sequences
        self.peptidecounter += len(pepkeys)

        for pep in pepkeys:
            # collect all the data for each peptide and create Peptide objects
            rx_pepdata = re.compile('^' + pep + '(?![0-9])')
            allpeptidedata = {}
            for x in peps4spec.keys():
                if rx_pepdata.search(x):
                    allpeptidedata[x] = peps4spec[x]

            # first create the peptide instance
            pepobj = peptide.ETpeptide(allpeptidedata, self)

            # then test if the sequence is valid
            if pepobj.isValidSequence():
                # only include peptides if the sequence is valid
                self.spectra[query]['numpeps'] += 1
                stats['numpeps'] += 1
                peplist.append(pepobj)
            else:
                stats['numfailedpeps'] += 1

        # now do the QC of the peptide set
        if peplist:
            self.doPeptideSetQC(peplist, 1)

        for pep in peplist:
            if pep.useinprot == 0 or pep.retain == 0:
                continue
            seq = pep.sequence
            score = pep.score

            # build dictionary of pep sequence to protein accession
            if seq in sequences:
                sequences[seq] += 1
                if pep.is_hook > seq2acc[seq]['hook']:
                    # replace non-hook data with hook data
                    seq2acc[seq]['hook'] = pep.is_hook
                    seq2acc[seq]['score'] = score
                    seq2acc[seq]['use'] = 1
                elif pep.is_hook == seq2acc[seq]['hook'] and score > seq2acc[seq]['score']:
                    # if hook status the same and score better use the new score
                    seq2acc[seq]['score'] = score

                if pep.czrank < seq2acc[seq]['bestczrank']:
                    seq2acc[seq]['bestczrank'] = pep.czrank
            else:
                sequences[seq] = 1
                seq2acc[seq] = dict(prots=pep.proteins[:], hook=pep.is_hook, score=score, bestczrank=pep.czrank)

    def doPeptideSetQC(self, peps, isETpep=0):
        '''
        @brief triggers the global processing of peptides and storage of global peptide parameters
        '''
        # referenced object variables
        analt = self.analtimes
        specs = self.spectra
        # peptides = self.peptides
        isModif = self.isModified
        scoreDelta = self.cfg.parameters['general']['delta']
        minHookLength = self.cfg.parameters['general']['minhooklength']

        qry = peps[0].query
        peps[0].czrank = 1
        seq = peps[0].sequence
        modsRelevant = peps[0].modsRelevant
        score = peps[0].score
        delta_mod = -1
        delta_seq = -1

        # calc modstrings, ranking and delta scores for the peptide set
        p = 1
        czrank = 1
        toprankseqs = [peps[0].sequence]
        peps[0].retain = 1
        lastscore = score
        while p < len(peps):
            if delta_seq == -1 and seq != peps[p].sequence:
                delta_seq = score - peps[p].score
            elif delta_mod == -1 and seq == peps[p].sequence and modsRelevant != peps[p].modsRelevant:
                delta_mod = score - peps[p].score

            if peps[p].score < lastscore:
                czrank += 1
                lastscore = peps[p].score
            peps[p].czrank = czrank

            if czrank == 1:
                peps[p].retain = 1
                if not peps[p].sequence in toprankseqs:
                    toprankseqs.append(peps[p].sequence)
            elif peps[p].sequence in toprankseqs:
                peps[p].retain = 1

            p += 1

        # copy da_delta scores to the spectrum dictionaries
        if not isModif(modsRelevant):
            # no mods so assign a score of -2000
            delta_mod = -2000
        elif delta_mod < 0:
            # only one sequence copy of that sequence
            delta_mod = score
        specs[qry]['delta_mod'] = delta_mod

        if delta_seq < 0:
            # only one sequence
            delta_seq = score
        specs[qry]['delta_seq'] = delta_seq

        if delta_seq > scoreDelta and len(toprankseqs[0]) >= minHookLength:
            p = 0
            # add hook to list for estimation of the real run time
            self.hookpeps.append(qry)
#            print '%s\t%f\t%i' % (peps[0].sequence,peps[0].mass,peps[0].query)

            for pep in peps:
                if pep.czrank == 1:
                    pep.is_hook = 1
                    self.hookppm.append((pep.da_delta, pep.da_delta / pep.mass * 1e6))
                else:
                    pep.useinprot = 0
        if isETpep:
            self.hdfMascot.writePeptides(peps, 'etpeptides')
        else:
            self.hdfMascot.writePeptides(peps, 'peptides')

        return

    def isModified(self, modstring):
        try:
            mod = int(modstring)
            if mod:
                return True
            else:
                return False
        except ValueError:
            return True

    def doPostPeptides(self):
        '''
        @brief triggers the global processing of peptides and storage of global peptide parameters
        '''

        # make sure the sequence to protein accession data is saved
        self.hdfMascot.writeSeq2Acc(self.seq2acc)

        stats = self.stats
        stats['numhooks'] = len(self.hookppm)

        self.logs.datlog.info('%d hookpeptides found' % stats['numhooks'])

        if stats['numhooks'] == 0:
            # no hook peptides
            stats['absdelta_avg'] = 0.0
            stats['reldelta_avg'] = 0.0
            stats['absdelta_stdev'] = 0.0
            stats['reldelta_stdev'] = 0.0
        else:
            ppm = np.array(self.hookppm)
            avg = np.average(ppm, 0)
            stdev = np.std(ppm, 0)
            stats['absdelta_avg'] = avg[0]
            stats['reldelta_avg'] = avg[1]
            stats['absdelta_stdev'] = stdev[0]
            stats['reldelta_stdev'] = stdev[1]

        self.hdfMascot.writeStatistics(stats)

    def findBestProtein(self):
        '''
        @brief use the seq2acc data to calculate the highest scoring protein
        '''

        protDict = {}
        bestScore = 0

        for seq in self.seq2acc:
            pepData = self.seq2acc[seq]
            protAccessions = [x['accession'] for x in pepData['prots']]
            score = pepData['hookscore']
            for acc in protAccessions:
                try:
                    protein = protDict[acc]
                    protein['score'] += score
                    protein['peps'] += pepData['numpep']
                    if pepData['hook']:
                        protein['hooks'] += 1
                except:

                    protDict[acc] = dict(score=score, peps=pepData['numpep'], hooks=0)
                    if pepData['hook']:
                        protDict[acc]['hooks'] = 1

                if protDict[acc]['score'] > bestScore:
                    bestScore = protDict[acc]['score']

        bestProts = []
        for acc in protDict:
            protein = protDict[acc]
            if protein['score'] == bestScore:
                bestProts.append(dict(accession=acc, hookScore=protein['score'],
                                      numHookPeps=protein['hooks'], numPeps=protein['peps']))

        bestProts.sort(key=lambda x: x['hookScore'])
        bestProt = bestProts[-1]
        # fetch peptide data and calculate peak capacity
        # fetch the peptide data
        peptides = self.hdfMascot.getPeptidesFromProtein(bestProt['accession'])

        # calculate peak capacity from best protein
        fwhmList = []
        rtList = []
        for seq in peptides:
            for id in peptides[seq]['peptides']:
                if id['rt'] > 0:
                    fwhmList.append(id['fwhm'])
                    rtList.append(id['rt'])
        stats = Statistics()
        try:
            capacity = 1 + (max(rtList) - min(rtList)) / stats.mean(fwhmList)
        except ValueError:
            capacity = 0
        bestProt['capacity'] = capacity
        bestProt['name'] = self.hdfMascot.getNameFromProteinAcc(bestProt['accession'])
        self.logs.datlog.info('Best protein: %s %s (hook score = %.1f)' % (bestProt['accession'], bestProt['name'],
                                                                           bestProt['hookScore']))
        self.hdfMascot.writeBestProtein(bestProt)
        return bestProt

    def addProteins(self, proteins):
        '''
        @brief processes the protein data for a single spectrum
        @param peps4spec <dictionary>: containing all peptides for one spectrum
        '''
        maxname = 0
        protList = []
        for acc in proteins:
            mass, name = proteins[acc].split(',', 1)
            protList.append(dict(accession=acc, mass=float(mass), name=name[1:-1]))
            proteins[acc] = dict(mw=float(mass), name=name[1:-1])
            if len(name) > maxname:
                maxname = len(name)

        if protList:
            self.hdfMascot.writeProteins(protList)

    def addQuerySpectra(self, query):
        '''
        @brief processes the query data including spectra
        @param query <dictionary>: containing the spectral data
        '''
        if query['query'] == '1':
            self.logs.datlog.info('Loading query spectral data')

        if 'title' in query:
            titledic = self.parseTitleData(query['title'])
        else:
            titledic = self.parseTitleData('')
        num = int(query['query'])

        if num in self.spectra:
            self.spectra[num].update(titledic)
            if 'msmsid'not in self.spectra[num]:
                self.logs.datlog.critical('no msmsid for spectrum number %s (%s) ' % (str(num), str(self.spectra[num])))
        else:
            self.spectra[num] = titledic
            if 'msmsid'not in self.spectra[num]:
                self.logs.datlog.critical('no msmsid for spectrum number %s (%s) ' % (str(num), str(self.spectra[num])))

    def doPostParsing(self):
        '''
        @brief do the post query parsing events
        '''
        # add check here
        self.logs.datlog.info('Writing query data')
        self.hdfMascot.writeQueries(self.spectra)
        self.hdfMascot.createTableIndexes(0)

    def doTimeAnalysis(self):
        '''
        @brief do the post query parsing events
        '''

        analt = self.analtimes
        specs = self.spectra

        # find last RT and group the early RTs
        for hook in self.hookpeps:
            rt = specs[hook]['rt']
            if rt > analt['last']:
                analt['last'] = rt
            if rt < analt['first']:
                analt['early'].append(rt)

        analt['early'].sort()
        ti = 0
        # find the first pair of RT's within 1 min as true start of elution
        while len(analt['early']) > 0 and analt['early'][0] - ti > 1:
            ti = analt['early'].pop(0)
        analt['first'] = ti
        analt['range'] = analt['last'] - analt['first']

    def addIndex(self, index):
        '''
        @brief writes index data to mascot.index table
        @param index <dictionary>: containing index data
        '''
        self.hdfMascot.writeIndex(index)

    def parseTitleData(self, title):
        '''
        @brief parses the title data from the mascot title string.  Mascot replaces punctuation
        '%2c' = ',', '%3a' = ':' and '%2e' = '.'
        @param title <string>: containing the spectrum TITLE data
        @return titles <dictionary>: containing the parsed values
        '''
        self.logs.datlog.debug('Parsing title')
        titles = {}
        expected = {'msmsid': '', 'rt': 0.0, 'surveyid': '', 'parent': 0.0, 'AnalTime': 0.0, 'Activation': None}

        if title == '':  # no title data so return minimum empty data set.
            return expected

        title = title.replace('%2e', '.')
        keydata = title.split('%2c')

        if len(keydata) == 1:  # for old mgf file (LCQ data)
            keydata = title.split('.')
            titles = expected.copy()

            titles['msmsid'] = 'F%0.6i' % keydata[2]

            return titles
        for data in keydata:
            pair = data.split('%3a')

            if len(pair) == 2:
                if pair[0] in ['msmsid', 'Activation']:
                    titles[pair[0]] = pair[1]
                elif pair[0] == 'survey':
                    titles['surveyid'] = pair[1]
                else:
                    try:
                        titles[pair[0]] = float(pair[1])
                    except:
                        self.logs.datlog.debug('cannot convert to float: %s' % str(pair[1]))
                        titles[pair[0]] = str(pair[1])
            elif len(pair) == 1:
                titles['msmsid'] = pair[0]
        self.logs.datlog.debug('titles: %s' % str(titles))

        # tidy up by making sure all the expected fields exist.
        for ex in expected:
            tmp = titles.setdefault(ex, expected[ex])

        return titles

    def parseSpectrumIons(self, ionstext):
        '''
        @brief parses the spectrum data by first splitting the string into ions then converting
        the mass-intensity strings into a tuple of floats (mass, intensity)
        @param ionstext <string>: with individual ions separated by , and mass-intensity values separated by :
        @return ionslist <list>: containing a list of the ion pairs sorted by mass
        '''

        ionslist = []
        templist = ionstext.split(',')
        self.logs.datlog.debug('Found %d ions' % len(templist))
        for ion in templist:
            data = ion.split(':')
            ionslist.append((float(data[0]), float(data[1])))
        ionslist.sort()
        return ionslist

    def doQuantitation(self):
        '''
        doQuantiation() -> void

        Reads the quantitation data for files that are iTRAQ or SILAC.
        '''

        # first check if the file is iTRAQ/SILAC
        rx_quant = re.compile('(_itraq|_silac){1}')
        file = self.parameterData['COM'].lower()
        match = rx_quant.search(file)
        if match:
            # build up a list of all the spectra to be quantified.
            self.logs.datlog.debug('File needs quantifying: %s' % self.parameterData['COM'])
            quandata = {}
            # specs = self.peptides.keys()
            specs = self.spectra.keys()
            specs.sort()
            for spec in specs:
                if self.spectra[spec]['origfile'] not in quandata:
                    quandata[self.spectra[spec]['origfile']] = {}
                quandata[self.spectra[spec]['origfile']][self.spectra[spec]['origspec']] = spec
            self.logs.datlog.debug('Quant Spectra: %s' % (str(quandata)))

    def doTopHitQC(self, hit):
        '''
        @brief performs QC on the top mascot protein
        @param hit <dictionary>: containing the protein data
        '''
        specs = self.spectra
        blank = (-1, 0, 0)
        if hit['score'] > 0:
            # do the protein QC
            print '  QC: %s %s (%d/%d)' % (hit['acc'], hit['name'], hit['score'], hit['peps'])
            peps = {}
            bsa = self.hdfMascot.getHitPeptides(hit['acc'])
            cover = {}
            data = np.ndarray(len(bsa), dtype=[('fwhm', float), ('ppm', float), ('rt', float)])
            for p in range(len(bsa)):
                seq = bsa[p]['sequence']
                specid = self.hdfMascot.getSpecidFromQuery(bsa[p]['query'])
                header = self.hdfMascot.getsMSMSHeaderFromSpecid(specid)

                drow = (header['fwhm'], bsa[p]['da_delta'] / bsa[p]['mass'] * 1000000, header['rtapex'])
                if bsa[p]['is_hook']:
                    data.put(p, drow)
                else:
                    data.put(p, blank)
                peps[seq] = max(peps.get(seq, 0), bsa[p]['score'])
                for pr in self.seq2acc[seq]['prots']:
                    if pr['accession'] == hit['acc']:
                        for i in range(pr['start'], pr['start'] + len(seq)):
                            cover[i] = 1
            data.sort(order='fwhm')
            i = 0
            while data[i]['fwhm'] == -1:
                i += 1
            data = data[i:]

            print '\tscore         = %d' % (np.sum(peps.values()))
            print '\tnum peptides  = %d' % (len(bsa))
            print '\tnum hook peps = %d' % len(data)
            print '\tcoverage      = %.1f %%' % (len(cover) / 6.07)
            print '\taccuracy (mDa)= %.3f +/- %.3f' % (np.mean(bsa['da_delta']) * 1000, np.std(bsa['da_delta']) * 1000)
            print '\taccuracy (ppm)= %.2f +/- %.2f' % (np.mean(data['ppm']), np.std(data['ppm']))
            print '\tRT range      = %.2f to %.2f min' % (self.analtimes['first'], self.analtimes['last'])
            print '\taverage width = %.1f s (%.1f - %.1f s)' % (np.mean(data['fwhm']) * 60, np.min(data['fwhm']) * 60,
                                                                np.max(data['fwhm']) * 60)
            print '\t   capacity   = %.1f' % (self.analtimes['range'] / np.mean(data['fwhm']))
            print '\tmedian width  = %.1f s' % (np.median(data['fwhm']) * 60)
            print '\t   capacity   = %.1f' % (self.analtimes['range'] / np.median(data['fwhm']))
