#!/usr/bin/env python

# python modules
from pathlib import Path
import os
import sys
import numpy as np

# CommonUtils imports

sys.path.insert(0, '..')
import CommonUtils.ExceptionHandler as ExHa
from CommonUtils.tools import Stopwatch
from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.LoggingManager import Logger
from CommonUtils.hdf5Results import HDF5Results
from CommonUtils.QuantMethodHandler import QuantMethods
import CommonUtils.progressbar as progBar
from CommonUtils.progressReport import progressReport

# application modules
from peptidesetmanager import PeptideSetManager
from hdfFile import HDFfile
from fileData import FileData


class DATimporter:
    def __init__(self, cfg, logger):
        """
        initialise the DATimporter
        @param cfg <ConfigManager object>: containing all the parameters from config file and comandline
        @return:
        """
        self.setsManager = PeptideSetManager(cfg, logger)
        self.hdfFiles = {}
        self.hdfFilesIdx = {}
        self.importData = FileData(None, None, logger.log)
        self.cfg = cfg
        self.pep2unique = {}

    @staticmethod
    def calculateFDR(FDRdata, threshold):
        scorekeys = sorted(FDRdata.keys(), reverse=True)
        cumrev = 0
        cumfwd = 0
        score2fdr = {}
        scoreatthreshold = 0
        for score in scorekeys:

            fwd, rev = FDRdata[score]
            logger.log.debug('fwd %s, rev %s @ score %s ' % (fwd, rev, score))
            cumfwd += fwd
            cumrev += rev
            try:
                fdr = cumrev / float(cumfwd)
            except ZeroDivisionError:
                logger.log.warning('Take a look at the data. There is something seriously wrong that you got no '
                                   'FWD hits')
                fdr = 1
            truehits = cumfwd - cumrev
            try:
                local_fdr = rev / float(fwd)
            except ZeroDivisionError:
                local_fdr = 1
            if fdr < threshold:
                scoreatthreshold = score
            logger.log.debug('fdr %s at score %s' % (fdr, score))

            score2fdr[score] = (cumfwd, cumrev, truehits, fdr,  local_fdr)
        logger.log.debug('finishing fdr calc with score at threshold %s' % scoreatthreshold)

        return scoreatthreshold, score2fdr

    def addHDFfile(self, hdfPath, idx, logger):
        """
        @brief adds an HDF5 file to the collection of files to process
        @param hdfPath <path>: path of the HDF5 file
        @return:
        """
        hdf = HDFfile(hdfPath, logger)

        # create fileData object to hold the relevant file summary data
        fd = FileData(hdfPath, hdf, logger.log)
        fd.indexname = '%s'.zfill(3) % idx
        self.hdfFiles[hdfPath.name] = fd
        self.hdfFiles[hdfPath.name].addMascotModifications(hdf.readModsTable())
        self.hdfFilesIdx[hdfPath.name] = 'f%s' % idx

        # look up the hook definition from the hdf5 file
        hdf.hdf.readOpen()
        data = hdf.hdf.getDataEqual('/%s/config' % hdf.importPath, 'parameter', 'delta')
        delta = int(data[0][2])
        data = hdf.hdf.getDataEqual('/%s/config' % hdf.importPath, 'parameter', 'minhooklength')
        minLength = int(data[0][2])
        decoyhitidentifier = hdf.hdf.getDataEqual('/%s/config' % hdf.importPath, 'parameter', 'decoyhitidentifier')
        decoysearchfromstart = hdf.hdf.getDataEqual('/%s/config' % hdf.importPath, 'parameter', 'decoysearchfromstart')
        decoyreplacementstring = hdf.hdf.getDataEqual('/%s/config' % hdf.importPath, 'parameter', 'decoyreplacementstring')

        # add information about how decoy hits are made identifiable

        try:
            self.setsManager.decoyhitidentifiers.add(decoyhitidentifier[0][2])
            self.setsManager.decoyreplacementstring = decoyreplacementstring[0][2]
            if decoysearchfromstart[0][2] == 'True' or decoysearchfromstart[0][2] == '1':
                self.setsManager.decoysearchfromstart = 1
        except:
            print 'problem trying to set decoyhitidentifer data: setting to default of ###REV###'
            self.setsManager.decoyhitidentifiers.add('###REV###')
            self.setsManager.decoysearchfromstart = 0
            self.setsManager.decoyreplacementstring = ''

        hdf.hdf.close()

        # read the new peptides from HDF5 and combine them with the existing data
        seq2acc = hdf.readSeq2accTable()
        proteins = hdf.readProteinTable()
        self.setsManager.proteinhitdata.update(proteins)
        self.setsManager.addSequenceData(seq2acc, proteins)
        allpeptides = hdf.readPeptides()
        self.setsManager.setFDRData(allpeptides)

    def finaliseProteins(self, peptidescoreatthreshold):
        """
        @brief controls the processing of the peptideSets: removal of peptideSets with too few hook peptides,
        then removal of duplicate peptideSets and optionally the removal of subsets.  All peptide in the Sets are
        tested for  their first use and the peptideSets have protein QC calculated.
        @return:
        """

        cfgGeneral = cfg.parameters['general']
        logger.log.info('%7i proteins from all peptides' % len(self.setsManager.peptideSets))
        logger.log.info('removing proteins with <%i hookpeptides' % cfgGeneral['minnumhook'])
        self.setsManager.finalisePeptideSets(cfgGeneral['minnumhook'], peptidescoreatthreshold)
        sw.rec('hook filter')

        logger.log.info('%7i proteins after hook peptide and unique filter' % len(self.setsManager.peptideSets))
        logger.log.info('removing duplicates')
        self.setsManager.removeDuplicates()

        sw.rec('remove duplicates')
        logger.log.info('%7i proteins after removal of duplicates' % len(self.setsManager.peptideSets))

        self.setsManager.calcIsFirstUseOfSequence()
        logger.log.info('%7i unique sequences flagged' % len(self.setsManager.peptide2set))
        sw.rec('is first use')
        logger.log.info('Calculating peptideSet QC values')

        self.setsManager.calcPeptideSetsQC(self.importData)
        if cfgGeneral['groupbygenename']:
            # will also group proteins by shared gene name
            self.setsManager.groupPerGene()

        # now calculate uniqueness using only FDR-threshold-passing peptides
        self.pep2unique = self.setsManager.calcUniqueness(peptidescoreatthreshold)
        self.setsManager.calcProteinFDR()
        sw.rec('peptide set QC')

    def updateHDF5(self):
        """
        @brief controls the updating of the data to the hdf5 results file
        @return finalMessage <string>: constructed from the protein data this is the RESULT stored in the DB
        """
        pep2unique = self.pep2unique
        baseContext = 'updateHDF5: '
        context = 'updateHDF5'
        try:
            # find the peptide sequences that are being imported
            usedPeps = self.setsManager.findUsedPeptides()
            logger.log.info('there are %s usedPeps' % len(usedPeps))

            context = baseContext + 'Retrieving sample IDs'

            sample_ids = range(1,  len(self.hdfFiles)+1)
            # create proteinset and proteinhit data
            starting_protein_group_no = 1
            self.setsManager.setProteinGroupNo(starting_protein_group_no)

            logger.log.info('adding protein group data to HDF5')

            logger.log.debug(str(self.hdfFiles.keys()))
            spectrum_id = 0
            peptide_id = 0
            hdfFileList = self.hdfFiles.keys()
            hdfFileList.sort()

            for key in hdfFileList:
                baseContext += '%s: ' % key
                logger.log.log(logger.PROCESS, 'Integrating Spectrum, Peptide & Quantification data from %s' % key)
                # collect fileData
                hdf = self.hdfFiles[key]
                hdfObj = hdf.hdfObject

                # set the current sample_id from the list of IDs extracted from the DB
                current_sample_id = sample_ids.pop()

                hdf.acquired_spectra, hdf.mascot_matched_spectra, numIsotopes, runTime = hdfObj.getNumbers()

                # read the Mascot data
                context = baseContext + 'Reading Mascot data'
                tmp = hdfObj.readImporterData(usedPeps, hdf)
                peptides = tmp[0]
                queryDict = tmp[1]
                headerArray = tmp[2]
                quantArray = tmp[3]
                quantExtraArray = tmp[4]

                # generate set of all spectra with quant
                quantSpectra = set()
                for q in quantArray:
                    quantSpectra.add(q['spec_id'])
                hdf.spectra_in_qc_proteins = len(peptides)

                logger.log.debug('getting spectrum_ids')
                context = baseContext + 'Retrieving spectrum IDs'

                acqTime, hdf.idAct, hdf.quanAct = hdfObj.getTimeAndActivation()
                # create blank lists to hold data for writing to hdf5 file
                spectrum_list = []
                peptide_list = []
                quant_list = []
                logger.log.info('collating spectrum, peptide & quant data')
                progRep = progressReport(len(queryDict), hdfObj.hdfFilePath.stem, 'collating query data', 'queries')
                for idx, q in enumerate(queryDict):
                    # loop round all the required spectra
                    progRep.next()
                    context = baseContext + 'query %i: Setting spectrum data' % q
                    # extract a spectrum_id from the list
                    spectrum_id += 1
                    query = queryDict[q]
                    spec = int(query['spec_id'])
                    context = baseContext + 'spectrum %i: Updating DB with spectrum data' % spec
                    # add spectrum data to spectrum_list
                    header = self.filterArrayEqual(headerArray, 'spec_id', spec)

                    spectrum_list.append(self.makeSpectrumDict(spectrum_id, current_sample_id, query, acqTime,
                                                               header))

                    # find the appropriate peptides
                    pepList = peptides[q]
                    logger.log.debug('there are %s in peplist %s' % (len(pepList), str(pepList)))
                    quantFound = 0

                    # this list will hold all peptides returned from makePeptideDictList and then filter
                    # those non-rank1 equivalents based on the score of the rank 1 peptide
                    tmplist = []

                    quantData = self.filterArrayEqual(quantArray, 'spec_id', spec)
                    if len(quantData) > 0:
                        hdf.quantified_spectra += 1
                    for pep in pepList:
                        # find the sets that the peptide belongs to and add to the peptide_list
                        sets = self.setsManager.peptide2set[pep['peptide']]
                        context = baseContext + 'spectrum %i: Creating peptide data entries for hdf5' % spec
                        tmp, qf = self.makePeptideDictList(spectrum_id, pep, query, sets, hdf, pep2unique, quantData)
                        tmplist.extend(tmp)
                        peptide_list += tmp
                        quantFound += qf

                    # only keep rank1 equivalent peptides (based on score)
                    tmplist.sort(key=lambda x: x['rank'])
                    toprankscore = tmplist[0]['score']
                    tmplist = [x for x in tmplist if x['score'] == toprankscore]

                    if quantMethID and quantFound:
                        # extract quantification data for the spectrum
                        if len(quantExtraArray) > 0:
                            quantExtraData = self.filterArrayEqual(quantExtraArray, 'spec_id', spec)
                        else:
                            quantExtraData = None
                        context = baseContext + 'spectrum %i: Creating quantification data entries for DB' % spec
                        newquant, deltas = self.gatherQuantData(spectrum_id, tmplist, header, quantData,
                                                                quantExtraData, quantSource)
                        quant_list += newquant

                        if quantSource == 'ms2':
                            context = baseContext + 'spectrum %i: Adding reporter ion delta data' % spec
                            hdf.addReporterDeltas(deltas)
                # pBar.finish()
                progRep.endReport()

                # calculate statistics
                context = baseContext + 'Calculating statistics'
                hdf.calcReporterStats()
                context = baseContext + 'Calculating delta m/z for fragment ions'

                context = baseContext + 'Updating sample table (%i)' % current_sample_id
                sample_data = hdf.getSampleDataDict(current_sample_id, key, runTime)

                hdf5results.writeSample(sample_data)

                self.importData.combineStatistics(hdf)

                # write data to HDF5
                context = baseContext + 'Updating spectrum table'
                logger.log.info('updating HDF5 with spectrum data')
                hdf5results.writeSpectrum(spectrum_list)

                if quantMethID:
                    context = baseContext + 'Updating specquant table'
                    logger.log.info('updating HDF5 with quant data')
                    hdf5results.writeSpecQuant(quant_list)

                context = baseContext + 'Retrieving peptide IDs'
                logger.log.info('updating HDF5 with peptide data')
                for pepdata in peptide_list:
                    pepdata['peptide_id'] = peptide_id
                    peptide_id += 1

                context = baseContext + 'Updating peptide table'
                hdf5results.writePeptide(peptide_list)
            hdf5results.createIndexes()

            logger.log.info('finalising HDF5 entries')
            hdf5results.writeFDRdata(self.importData.score2fdr, 'peptide')
            hdf5results.writeFDRdata(self.importData.proteinscore2fdr, 'protein')

            topScoringProteinInfo = self.setsManager.addPeptideSetDBdata(hdf5results, self.importData.proteinscore2fdr)
            runtimedata = self.importData.getSummaryStatisticsDict()

            hdf5results.writeStatistics(runtimedata)

            finalMessage = 'queries matched: %i / %s (%.1f%%) ' % (runtimedata['spectra_in_qc_proteins'],
                                                                   runtimedata['mascot_matched_spectra'],
                                                                   (runtimedata['spectra_in_qc_proteins'] /
                                                                    float(runtimedata['mascot_matched_spectra'])) * 100)
            finalMessage += 'spectra quantified: %i top hit %s (%s) ' % (runtimedata['quantified_spectra'], '', '')
            finalMessage += 'with total score %f and %i matched peptides (hook AND non hook)' % \
                            (topScoringProteinInfo[0], topScoringProteinInfo[2])

            baseContext = 'updateHDF5: '
            context = baseContext + 'Finalising HDF5 entries'
        except Exception, genEx:
            # make sure that there aren't any permanent changes
            ExHa.addContext(genEx, context)
            finalMessage = 'Error: %s' % ExHa.oneLineRepr(genEx)
            raise

        return finalMessage

    @staticmethod
    def filterArrayEqual(array, column, value):
        if len(array) > 0:
            arrayfilter = (array[column] == value)
            return array[arrayfilter]
        else:
            return np.array([])

    @staticmethod
    def makeSpectrumDict(spectrum_id, sample_id, query, acqTime, header):

        specDict = dict(spectrum_id=spectrum_id,
                        msms_id=query['msms_id'], query=int(query['query']),
                        neutral_mass=query['neutral_mass'], aquisitiontime=acqTime, charge_state=header[0]['charge'],
                        precursor_mz=header[0]['precmz'], start_time=header[0]['rt'],
                        survey_id=int(header[0]['survey_spec']), parent_ion=header[0]['setmass'], rt=header[0]['rt'],
                        peak_rt=header[0]['rtapex'], s2i=header[0]['s2i'], peak_fwhm=header[0]['fwhm'],
                        fragenergy=header[0]['fragenergy'], precursor_noise=header[0]['thresh'],
                        sumreporterarea=header[0]['sumreparea'], sumreporterint=header[0]['sumrepint'],
                        peak_inten_raw=header[0]['inten'], peak_inten_smoothed=header[0]['smooth'],
                        peak_inten_integrated=header[0]['integ'], peak_area=header[0]['area'],
                        precursor_area=header[0]['c12'])

        specDict['sample_id'] = sample_id

        return specDict

    def makePeptideDictList(self, spectrum_id, pep, query, sets, hdf, pep2unique, quantData):

        score2fdr = self.importData.score2fdr
        try:
            is_unique = pep2unique[pep['peptide']]
        except KeyError:
            is_unique = 0
        pepDict = dict(spectrum_id=spectrum_id, score=pep['score'],
                       rank=pep['rank'], peptide=pep['peptide'],
                       is_unique=is_unique, failed_fdr_filter=pep['failed_fdr_filter'],
                       mw=pep['mw'], da_delta=pep['da_delta'],
                       is_hook=pep['is_hook'], is_quantified=0,
                       is_decoy=0, missed_cleavage_sites=pep['missed_cleavage_sites'],
                       delta_seq=query['delta_seq'], ppm_error=pep['da_delta'] / pep['mw'] * 1e6,
                       delta_mod=query['delta_mod'], in_protein_inference=0,
                       fdr_at_score=score2fdr[round(pep['score'])][3])
        mods = hdf.modsHandler.generateModifications(pep)
        quantFound = 0
        # here we check that the modification relating to the quantification label is in the
        # mascot-assigned modifications
        if quantMethID and len(quantData) > 0:
            if isNullQuant or quantMethData['mascot_name'] in mods['positional_modstring']:
                quantFound = 1
                pepDict['is_quantified'] = 1
        elif quantData and quantSource == 'ms1':
            quantFound = 1
            pepDict['is_quantified'] = 1
        logger.log.debug('isquantified is %s for peptide %s (len quantData %s) ' %
                         (quantFound, pep['peptide'], len(quantData)))
        pepDict.update(mods)
        pepDictList = []
        for acc in sets:
            pepSet = self.setsManager.peptideSets[acc]
            seq_start = pepSet.peptideData[pep['peptide']]['seq_start']
            seq_end = pepSet.peptideData[pep['peptide']]['seq_end']

            if pepSet.is_decoy:
                pepDict['is_decoy'] = 1
            if pepDict['peptide'] in pepSet.validpeptides:
                pepDict['in_protein_inference'] = 1
            pepDict['seq_start'] = seq_start
            pepDict['seq_end'] = seq_end
            pepSet.adddataFromSet(pepDict)
            pepDictList.append(pepDict.copy())

        return pepDictList, quantFound

    def gatherQuantData(self, spectrum_id, pepData, headerref, quantData, quantExtraData, quantSource):
        """
        This method gathers and prepares quant data from different sources. These data will later be used to filter
        and select quant values based on user-supplied parameters
        :param spectrum_id: internally generated spectrum id to which quant values will be linked
        :param pepData: data relating to peptide sequence and associated information
        :param headerref: data from msmsheader table (supplying data relating to original MS/MS event)
        :param quantData: numpy array with quant data for spectrum id
        :param quantExtraData: either numpy array or None
        :param quantSource: which source quant data are from (MS1 or MS2)
        :return: quantList list of dictionaries of concatenated quant-values and filterable values.
        List is only more than one element when quant values are in more than one protein group, deltas keeps the
        reporter ion mz deltas
        """
        score2fdr = self.importData.score2fdr

        pepData.sort(key=lambda x: x['rank'])
        if len([x for x in pepData if x['rank'] == 1]) > 1:
            is_unique = 0
        else:
            # take lowest unique value as an approximation for true uniqueness.
            is_unique = min(p['is_unique'] for p in pepData)
        firstpep = pepData[0]
        logger.log.debug('started gatherQuantData for %s, len(quantData) %s' % (firstpep, len(quantData)))
        score = round(firstpep['score'])
        fdr_at_score = score2fdr[score][3]
        delta_seq = firstpep['delta_seq']
        peptide_length = len(firstpep['peptide'])
        quantList = []

        protein_group_nos = set([x['protein_group_no'] for x in pepData])
        deltas = []
        if len(quantData) > 0:
            # accumulate numbers
            s2i = headerref['s2i']
            p2t = headerref['c12'] / headerref['thresh']

            if len(protein_group_nos) > 1:
                logger.log.debug('more than one proteinset (%s) !! %s' % (str(protein_group_nos), spectrum_id))
            # all data in quantData should be based on the rank 1 peptide sequence. If there are more than one entry
            # it's either due to membership of  more than one proteinset or peptide of different rank but the same
            # score and sequence composition.
            for q in quantData:
                for protein_group_no in protein_group_nos:
                    least_squares = 0
                    prior_ion_ratio = 0
                    source = ''
                    if 'ms1' in quantSource:
                        quant_raw = q['inten']
                        if quantExtraData is not None:
                            least_squares = max(quantExtraData['max_least_squares'])
                            prior_ion_ratio = max(quantExtraData['prior_ion_ratio_final'])
                            source = quantExtraData['source'][0]
                    else:
                        quant_raw = q['area']
                    specQuantDict = dict(protein_group_no=protein_group_no, spectrum_id=spectrum_id,
                                         delta_seq=delta_seq, is_unique=is_unique,
                                         isotopelabel_id=int(q['isolabel_id']), score=firstpep['score'],
                                         quant_raw=quant_raw, quant_isocorrected=q['corrected'],
                                         quant_allcorrected=q['corrected'], p2t=p2t,
                                         peptide_length=peptide_length, s2i=s2i,
                                         in_quantification_of_protein=0, fdr_at_score=fdr_at_score,
                                         least_squares=least_squares, prior_ion_ratio=prior_ion_ratio,
                                         ms1source=source)
                    quantList.append(specQuantDict.copy())
                deltas.append(q['mzdiff'])
        return quantList, deltas


if __name__ == '__main__':
    logger = None
    cfg = None

    installdir = Path(os.path.dirname(__file__))
    confile = installdir / Path('proteininference.cfg')
    cfg = ConfigManager(str(confile))
    ret = cfg.evaluateCommandLineArgs(sys.argv)
    searches = sorted(cfg.parameters['runtime']['filelist'])
    dataDir = cfg.parameters['runtime']['datadir']
    resultfile = dataDir.joinpath(cfg.parameters['runtime']['resultfilename'])

    try:
        # start the logging process
        logParam = cfg.parameters['logging']
        logPath = Path(dataDir.joinpath(logParam['logdir']))
        if not logPath.exists():
            logPath.mkdir(parents=True)
        logFile = logPath.joinpath(logParam['logfile'])

        logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], True)
        logger.setProtInferenceLogs()

        if resultfile.exists():
            resultfile.unlink()
        expID = None

        quantMethID = cfg.parameters['runtime']['quantmethod_id']
        isNullQuant = False
        if quantMethID:
            quantMethData = QuantMethods().getMethodByID(quantMethID)
            quantSource = quantMethData['source']
            for id in quantMethData['quantmasses']:
                for isotope in quantMethData['quantmasses'][id]:
                    if isotope['mass'] == 0:
                        isNullQuant = True
        else:
            quantMethData = {}
            quantSource = ''

        hdf5results = HDF5Results(str(resultfile))
        hdf5results.appendOpen()
        hdf5results.createTables()
        hdf5results.addConfigParameters(cfg.parameters, 'postMascot', '2 proteinInference')

        archives = []
        logger.log.info('Results file = %s' % resultfile.name)
        logger.log.info('MinNumHooks  = %i' % cfg.parameters['general']['minnumhook'])

        sw = Stopwatch()

        allSearchData = {}
        fractions = ''

        if len(searches) > 1:
            isMerged = True
            msg = 'merged sample'
        else:
            isMerged = False
            msg = 'single sample'

        importer = DATimporter(cfg, logger)
        logger.log.info(msg)
        searchDict = {}
        for idx, hdf5file in enumerate(searches):
            searchData = dict()
            searchData['hdf5name'] = hdf5file
            searchData['archivepath'] = cfg.parameters['runtime']['datadir']
            fullPath = searchData['archivepath']
            fi = fullPath.joinpath(searchData['hdf5name'])
            if not fi.exists():
                raise ExHa.FileNotFoundException('Missing file: %s' % str(fi))
            logger.log.info('Reading: %s' % str(fi))
            importer.addHDFfile(fi, idx+1, logger)
        sw.rec('file loading')
        fdrthreshold = cfg.parameters['general']['fdrthreshold']
        peptidescoreatthreshold, score2fdr = importer.calculateFDR(importer.setsManager.FDRdata, fdrthreshold)
        importer.importData.score2fdr = score2fdr

        importer.finaliseProteins(peptidescoreatthreshold)
        proteinscoreatthreshold, score2fdr = importer.calculateFDR(importer.setsManager.proteinFDRdata, fdrthreshold)

        importer.importData.proteinscore2fdr = score2fdr

        finalMessage = importer.updateHDF5()
        finalMessage = 'finished: %s' % resultfile.name

        sw.rec('update DB')

        logger.log.info('Total Spectra         = %6i' % importer.importData.acquired_spectra)
        logger.log.info('Total Spectra w Pep   = %6i' % importer.importData.mascot_matched_spectra)
        logger.log.info('Total Spectra Matched = %6i' % importer.importData.spectra_in_qc_proteins)
        logger.log.info('Total Spectra w Quant = %6i' % importer.importData.quantified_spectra)
        logger.log.info('Total Spectra All Rep = %6i' % importer.importData.numSpectraAllReporters)

        times = sw.stop()
        logger.log.info(sw.format())

        hdf5results.close()

    except ExHa.UsageError as useEx:
        ExHa.reformatException(useEx)
        logger.log.warning(ExHa.oneLineRepr(useEx))
    except Exception as genEx:
        # error
        if logger:
            logger.log.warning(ExHa.oneLineRepr(genEx))
            print ExHa.multiLineRepr(genEx)
        else:
            print ExHa.multiLineRepr(genEx)
        if cfg:
            ExHa.exportError2File(genEx, dataDir.joinpath(resultfile.stem + '.error'))
        else:
            ExHa.exportError2File(genEx, dataDir.joinpath('errors.error'))

        sys.exit(ExHa.oneLineRepr(genEx))
