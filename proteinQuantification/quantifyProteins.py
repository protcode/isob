# !/usr/bin/env python

"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This is a set of methods to perform protein fold change quantification based
reporter ion signals filtered according to certain criteria

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here:
https://github.com/protcode/isob/
"""


# python imports
from pathlib import Path
import sys
import numpy as np
import numpy.random as nr
sys.path.insert(0, '..')

# CommonUtils imports
from CommonUtils.ConfigManager import ConfigManager
import CommonUtils.ExceptionHandler as ExHa
from CommonUtils.LoggingManager import Logger
from CommonUtils.hdf5Results import HDF5Results
import CommonUtils.progressbar as progBar


class QuantifyProtController:
    def __init__(self,  hdf5quantprot, ncqcfg):
        self.hdf5File = hdf5quantprot.hdfFilePath
        self.hdf5quantprot = hdf5quantprot
        self.cfg = ncqcfg

    def getFilteredSpectrumids(self):
        """
        @brief get all spectum ids which pass filters in config file
        We basically create a intersection of spectra ids which pass the thresholds and return this
        @return set validspectrumids
        """
        p2tfilter = 'p2t > %s' % self.cfg.parameters['general']['p2tthreshold']
        s2ifilter = 's2i > %s' % self.cfg.parameters['general']['s2ithreshold']
        mascotfilter = 'score > %s' % self.cfg.parameters['general']['mascotthreshold']
        fdrfilter = 'fdr_at_score < %s' % self.cfg.parameters['general']['fdrthreshold']
        uniquefilter = 'is_unique==1'
        peplengthfilter = 'peptide_length >= %s' % self.cfg.parameters['general']['peplengthfilter']
        deltaseqfilter = '(delta_seq > %s) | (delta_seq == 0)' % self.cfg.parameters['general']['deltaseqfilter']
        validspectrumids = self.hdf5quantprot.getFilterSpectra(p2tfilter)
        self.cfg.log.debug('there are %s ids after p2tfilter' % len(validspectrumids))
        validspectrumids = validspectrumids & self.hdf5quantprot.getFilterSpectra(s2ifilter)
        self.cfg.log.debug('there are %s ids after p2tfilter & s2ifilter ' % len(validspectrumids))
        validspectrumids = validspectrumids & self.hdf5quantprot.getFilterSpectra(mascotfilter)
        self.cfg.log.debug('there are %s ids after p2tfilter & s2ifilter & mascotfilter ' % len(validspectrumids))
        validspectrumids = validspectrumids & self.hdf5quantprot.getFilterSpectra(peplengthfilter)
        validspectrumids = validspectrumids & self.hdf5quantprot.getFilterSpectra(uniquefilter)
        self.cfg.log.debug('there are %s ids after p2tfilter & s2ifilter & mascotfilter & peplengthfilter' %
                           len(validspectrumids))
        validspectrumids = validspectrumids & self.hdf5quantprot.getFilterSpectra(deltaseqfilter)
        self.cfg.log.debug('there are %s ids after p2tfilter & s2ifilter & mascotfilter & deltaseqfilter '
                           '& peplengthfilter ' % len(validspectrumids))
        validspectrumids = validspectrumids & self.hdf5quantprot.getFilterSpectra(fdrfilter)
        self.cfg.log.debug('there are %s ids after p2tfilter & s2ifilter & mascotfilter & deltaseqfilter & fdrfilter '
                           '& peplengthfilter ' % len(validspectrumids))

        return validspectrumids

    def getValidProteinGroups(self):
        """
        @brief Get all proteingroup ids from results table
        @return set proteingroup nos (identifiers for protein group)
        """
        specquantable = self.hdf5quantprot.getMS2Data()
        protein_group_nos = set([x['protein_group_no'] for x in specquantable])
        return protein_group_nos

    def performSimpleSumQuant(self, protein_group_nos, reference):
        """
        @brief perform simple sum ratio calculation using valid quant data from all spectra
        """
        all_proteins_quantdata = self.hdf5quantprot.getAllProteinDatafromSpecQuant()
        peptidedatafromsets = self.hdf5quantprot.getPeptideDataforSets()
        protein2quantdata = {}  # keeps fold change result,sum ion area, isotopelabel
        usedIsotopes = sorted(list(set([x['isotopelabel_id'] for x in all_proteins_quantdata])))
        pBar = progBar.ProgressBar(widgets=progBar.name_widgets, maxval=len(protein_group_nos),
                                   name='bootstrap quant').start()

        # scan through each protein group and perform quantification (according to method given) on quant values
        #
        for idx, protein_group_no in enumerate(protein_group_nos):
            pBar.update(idx)
            datadict = {}
            self.cfg.log.debug('starting for proteingroup %s' % protein_group_no)
            protlocation = all_proteins_quantdata['protein_group_no'] == protein_group_no
            data = all_proteins_quantdata[protlocation]
            peptidedataforset = peptidedatafromsets[protein_group_no]

            quantuniquepeps = set([peptidedataforset[spectrum_id] for spectrum_id in data['spectrum_id']])

            refdata = data[data['isotopelabel_id'] == reference]['quant_allcorrected']
            sumrefdata = refdata.sum()
            qupm = len(quantuniquepeps)
            refspectraquantified = len(refdata)
            if sumrefdata:
                datadict[reference] = [sumrefdata, refspectraquantified, qupm, (1, 0, 0)]
            else:
                datadict[reference] = [sumrefdata, refspectraquantified, qupm, (-1, -1, -1)]
                for isotopelabel_id in usedIsotopes:
                    datadict[isotopelabel_id] = [0, 0, qupm, (-1, -1, -1)]
                    # cannot continue as there are no references!
                    break
            for isotopelabel_id in usedIsotopes:

                self.cfg.log.debug('assessing isotopelabel %s' % isotopelabel_id)
                if isotopelabel_id != reference:
                    queryvalues = data[data['isotopelabel_id'] == isotopelabel_id]['quant_allcorrected']
                    sumqueryvalues = queryvalues.sum()
                    qssm = len(queryvalues)

                    self.cfg.log.debug('length ref %s, length query %s' % (len(refdata), len(queryvalues)))
                    if queryvalues.any() and refdata.any():
                        ratio_result = sumqueryvalues / sumrefdata
                        result = (ratio_result, -1, -1)
                        self.cfg.log.debug('simple ratio result for %s : %s ' % (protein_group_no, result))
                    elif refdata.any():
                        self.cfg.log.debug('there are no valid query values for %s: FC will be zero, as reference'
                                           ' is present' % protein_group_no)
                        result = (0, -1, -1)
                    else:
                        result = (-1, -1, -1)
                    datadict[isotopelabel_id] = [sumqueryvalues, qssm, qupm, result]
            protein2quantdata[protein_group_no] = datadict
        pBar.finish()
        return protein2quantdata

    def performBootstrapQuant(self, protein_group_nos,  reference):
        """
        @brief get all filtered spectra from spec quant table then fetch peptide data
        (sequence and spectrum id ) from all protein sets.
        for every protein group perform fold change calculation using bootstrap method
        @param protein_group_nos list of protein group ids
        @param reference id of value used for fold change calculation using bootstrap model
        """

        all_proteins_quantdata = self.hdf5quantprot.getAllProteinDatafromSpecQuant()
        peptidedatafromsets = self.hdf5quantprot.getPeptideDataforSets()
        protein2quantdata = {}  # keeps fold change result,sum ion area, isotopelabel
        usedIsotopes = sorted(list(set([x['isotopelabel_id'] for x in all_proteins_quantdata])))
        pBar = progBar.ProgressBar(widgets=progBar.name_widgets, maxval=len(protein_group_nos),
                                   name='bootstrap quant').start()
        for idx, protein_group_no in enumerate(protein_group_nos):
            pBar.update(idx)
            missingrefevents = 0
            datadict = {}
            self.cfg.log.debug('starting for proteingroup %s' % protein_group_no)
            protlocation = all_proteins_quantdata['protein_group_no'] == protein_group_no
            data = all_proteins_quantdata[protlocation]
            peptidedataforset = peptidedatafromsets[protein_group_no]

            quantuniquepeps = set([peptidedataforset[spectrum_id] for spectrum_id in data['spectrum_id']])

            totalquantevents = len(set(data['spectrum_id']))
            refdata = data[data['isotopelabel_id'] == reference]['quant_allcorrected']
            sumrefdata = refdata.sum()
            refspectraquantified = len(refdata)
            if refspectraquantified != totalquantevents:
                missingrefevents = totalquantevents - refspectraquantified
            qupm = len(quantuniquepeps)
            if sumrefdata:
                datadict[reference] = [sumrefdata, refspectraquantified, qupm, (1, 0, 0)]
            else:
                # if there are no quantified peptides from the reference label then we cannot calculate a fold change
                datadict[reference] = [sumrefdata, refspectraquantified, qupm, (-1, -1, -1)]
            for isotopelabel_id in usedIsotopes:
                missingqueryevents = 0
                self.cfg.log.debug('assessing isotopelabel %s' % isotopelabel_id)
                if isotopelabel_id != reference:
                    queryvalues = data[data['isotopelabel_id'] == isotopelabel_id]['quant_allcorrected']
                    sumqueryvalues = queryvalues.sum()
                    qssm = len(queryvalues)

                    if qssm != totalquantevents:
                        missingqueryevents = totalquantevents - qssm
                    self.cfg.log.debug('length ref %s, length query %s' % (len(refdata), len(queryvalues)))
                    if queryvalues.any() and refdata.any():
                        minquantspectra = cfg.parameters['general']['minquantspectra']
                        result = self.makeBootstrap(queryvalues.tolist()+[0]*missingqueryevents,
                                                    refdata.tolist()+[0]*missingrefevents, minquantspectra)
                        self.cfg.log.debug('bootstrap result for %s : %s ' % (protein_group_no, result))
                    elif refdata.any():
                        self.cfg.log.debug('there are no valid query values for %s: FC will be zero, as reference'
                                           ' is present' % protein_group_no)
                        result = (0, -1, -1)
                    else:
                        result = (-1, -1, -1)
                    datadict[isotopelabel_id] = [sumqueryvalues, qssm, qupm, result]
            protein2quantdata[protein_group_no] = datadict
        pBar.finish()
        return protein2quantdata

    def updateValidSpectra(self, validspectra):
        """
        @brief set in_quantification_of_protein in specquant table to 1 for  all spectra in .hdf5 file
        @validspectra set all validated spectrum ids
        """
        self.hdf5quantprot.updatein_quantification_of_protein(validspectra)

    def makeBootstrap(self, queryvalues, refvalues, minquantspectra):
        """
        @brief actually perform boot strap for isotope label for one protein provided the minimum spectra number
        are present
        @queryvalues list of all query reporter signals
        @refvalues list of all refvalues reporter signals
        @minquantspectra int minimum count of  non-null reporter signals for bootstrap calculation to proceed.
        """
        if len([x for x in queryvalues if x > 0]) < minquantspectra \
                and len([x for x in refvalues if x > 0]) < minquantspectra:
            self.cfg.log.debug('will not bootstrap for under %s valid spectra' % minquantspectra)
            try:
                # with satatement suppresses possible warnings about divide by zero or invalid numbers
                with np.errstate(divide='ignore', invalid='ignore'):
                    sumratio = sum(queryvalues) / sum(refvalues)
            except ZeroDivisionError:
                self.cfg.log.info('cannot perform sum ratio calcultion either')
                return -1, -1, -1
            return sumratio, -1, -1
        self.cfg.log.debug('started makeBootstrap with %s and %s' % (str(queryvalues), str(refvalues)))
        valcount = min(len(queryvalues), len(refvalues))
        nbs = 5000  # number of bootstrap iterations
        samples = nr.randint(0, valcount, (valcount, nbs))
        queryarray = np.array(queryvalues)[samples]
        refarray = np.array(refvalues)[samples]
        refsum = np.sum(refarray, 0)
        testsum = np.sum(queryarray, 0)
        with np.errstate(divide='ignore', invalid='ignore'):
            bootstrap = testsum / refsum
        t = np.isinf(bootstrap)
        bootstrap[t] = max(bootstrap[t == 0])
        result = self.quantile(bootstrap, [0.50, 0.025, 0.975])
        self.cfg.log.debug('result of bootstrap %s' % str(result))
        return result.tolist()

    @staticmethod
    def quantile(x, probs):
        """
        @brief compute quantiles defined in probs (quantiles in [0,1]) on an array x
        """
        if isinstance(x, list):
            x = np.array(x)                           # make sure we are dealing with an array
        n = len(x)

        probs = np.array(probs)                       # convert probs to array -> allows indexing

        index = (n - 1) * probs
        # with 5000 index is always: [2499.5 (50%), 124.975 (first 2.5%), 4874.025 (last 2.5% (97.5%))]
        lo = np.floor(index)
        lo = lo.astype('int')                         # conversion required to use as index
        # for 5000 lo is always [2499,  124, 4874]
        hi = np.ceil(index)
        hi = hi.astype('int')                         # conversion required to use as index
        # for 5000 hi is always [2500  125 4875]
        x.sort()

        i = index > lo                                # identify which indexes are problematic

        qs = x[lo]                                    # get quantiles

        h = (index - lo)[i]                           # compute gamma correction

        if len(h) > 0:                                # in case some
            qs[i] = (1 - h) * qs[i] + h * x[hi[i]]    # apply correction for problematic index

        return qs

    def saveQuantData(self, allfoldchangedata, reference):
        prepareddata = []
        for protein_group_no, data in allfoldchangedata.iteritems():
            for isotopelabel_id, isospecificdata in data.iteritems():
                sumdata = isospecificdata[0]
                qssm = isospecificdata[1]
                qupm = isospecificdata[2]
                protein_fold_change = isospecificdata[3][0]
                lower_confidence_level = isospecificdata[3][1]
                upper_confidence_level = isospecificdata[3][2]
                prepareddata.append(dict(protein_group_no=protein_group_no,
                                         isotopelabel_id=isotopelabel_id,
                                         reference_label=str(reference),
                                         protein_fold_change=protein_fold_change,
                                         lower_confidence_level=lower_confidence_level,
                                         upper_confidence_level=upper_confidence_level,
                                         sum_quant_signal=sumdata,
                                         qssm=qssm,
                                         qupm=qupm))

        self.hdf5quantprot.writeProteinQuant(prepareddata)

if __name__ == "__main__":
    logger = 0
    cfg = ConfigManager('./proteinquantification.cfg')
    try:
        configErr = cfg.evaluateCommandLineArgs(sys.argv)
        dataDir = cfg.parameters['runtime']['datadir']
        # start the logging process
        logParam = cfg.parameters['logging']
        logPath = Path(dataDir.joinpath(logParam['logdir']))
        if not logPath.exists():
            logPath.mkdir(parents=True)
        logFile = logPath.joinpath(logParam['logfile'])

        logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], False)
        logger.setLog('proteinQuantification')

        cfg.log = logger.log

        hdf5files = [Path(cfg.parameters['runtime']['datadir'], cfg.parameters['runtime']['filename'])]
        counter = 0
        # go through all results.hdf5 files in list: normally we should only have one for merged datasets
        for hdf5Path in hdf5files:
            logger.log.info('hdf5Path %s' % hdf5Path)
            logger.log.info('starting')
            hdf5quanprot = HDF5Results(hdf5Path)
            hdf5quanprot.appendOpen()
            hdf5quanprot.createQuantTables()
            hdf5quanprot.samplename = hdf5Path.name
            hdf5quanprot.cfg = cfg
            hdf5quanprot.addConfigParameters(cfg.parameters, 'postMascot', '5 quantifyProteins')
            QuantProt = QuantifyProtController(hdf5quanprot, cfg)
            counter += 1
            validspectra = QuantProt.getFilteredSpectrumids()
            QuantProt.updateValidSpectra(validspectra)
            protein_group_nos = QuantProt.getValidProteinGroups()
            referencelabel = cfg.parameters['runtime']['reference']
            if cfg.parameters['general']['quantmethod'] == 'bootstrap':
                proteingroup2foldchanges = QuantProt.performBootstrapQuant(protein_group_nos,
                                                                           cfg.parameters['runtime']['reference'])
            elif cfg.parameters['general']['quantmethod'] == 'simpleratio':
                proteingroup2foldchanges = QuantProt.performSimpleSumQuant(protein_group_nos,
                                                                           cfg.parameters['runtime']['reference'])
            else:
                message = 'quantmethod %s not yet supported' % cfg.parameters['general']['quantmethod']
                sys.exit(message)

            QuantProt.saveQuantData(proteingroup2foldchanges, cfg.parameters['runtime']['reference'])
            hdf5quanprot.close()

    except ExHa.czException as czEx:
        ExHa.reformatException(czEx)
        ExHa.addContext(czEx, 'Error during QuantifyProteins run')
        ExHa.exportError2File(czEx, cfg.parameters['runtime']['datadir'] / Path('errors.error'))
        if logger:
            logger.log.warning(ExHa.oneLineRepr(czEx))
        else:
            print ExHa.multiLineRepr(czEx)

    except Exception as genEx:

        ExHa.reformatException(genEx)
        ExHa.addContext(genEx, 'Error during QuantifyProteins run')
        ExHa.exportError2File(genEx, cfg.parameters['runtime']['datadir'] / 'errors.error')
        if logger:
            logger.log.warning(ExHa.oneLineRepr(genEx))
        else:
            print ExHa.multiLineRepr(genEx)
