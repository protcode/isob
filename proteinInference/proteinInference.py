#!/usr/bin/env python

# python modules
from pathlib import Path
import os
import sys

# CommonUtils imports

sys.path.insert(0, '..')
import CommonUtils.ExceptionHandler as ExHa
from CommonUtils.tools import Stopwatch
from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.LoggingManager import Logger
from CommonUtils.hdf5Results import HDF5Results
from CommonUtils.QuantMethodHandler import QuantMethods
import CommonUtils.progressbar as progBar

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
        hdf.hdf.close()

        # read the new peptides from HDF5 and combine them with the existing data
        seq2acc = hdf.readSeq2accTable()
        self.setsManager.proteinhitdata.update(hdf.readProteinTable())
        self.setsManager.addSequenceData(seq2acc)
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
                quanArray = tmp[3]

                hdf.spectra_in_qc_proteins = len(peptides)

                logger.log.debug('getting spectrum_ids')
                context = baseContext + 'Retrieving spectrum IDs'

                acqTime, hdf.idAct, hdf.quanAct = hdfObj.getTimeAndActivation()
                # create blank lists to hold data for writing to hdf5 file
                spectrum_list = []
                peptide_list = []
                quant_list = []
                logger.log.info('collating spectrum, peptide & quant data')
                pBar = progBar.ProgressBar(widgets=progBar.name_widgets, maxval=len(queryDict),
                                           name='collate data').start()
                for idx, q in enumerate(queryDict):
                    # loop round all the required spectra
                    pBar.nextPrimary()
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
                    for pep in pepList:
                        # find the sets that the peptide belongs to and add to the peptide_list
                        sets = self.setsManager.peptide2set[pep['peptide']]
                        context = baseContext + 'spectrum %i: Creating peptide data entries for hdf5' % spec
                        tmp, qf = self.makePeptideDictList(spectrum_id, pep, query, sets, hdf, pep2unique)
                        tmplist.extend(tmp)
                        peptide_list += tmp
                        quantFound += qf

                    # only keep rank1 equivalent peptides (based on score)
                    tmplist.sort(key=lambda x: x['rank'])
                    toprankscore = tmplist[0]['score']
                    tmplist = [x for x in tmplist if x['score'] == toprankscore]

                    if quantMethID and quantFound:
                        # extract quantification data for the spectrum
                        context = baseContext + 'spectrum %i: Creating quantitation data entries for DB' % spec
                        newquant, deltas = self.makeQuantDictLists(spectrum_id, spec,  tmplist, header, quanArray, hdf)

                        quant_list += newquant

                        if quantSource == 'ms2':
                            context = baseContext + 'spectrum %i: Adding reporter ion delta data' % spec
                            hdf.addReporterDeltas(deltas)
                pBar.finish()

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
        arrayfilter = (array[column] == value)
        return array[arrayfilter]

    @staticmethod
    def makeSpectrumDict(spectrum_id, sample_id, query, acqTime, header):

        specDict = dict(spectrum_id=spectrum_id,
                        msms_id=query['msms_id'], query=int(query['query']),
                        neutral_mass=query['neutral_mass'], aquisitiontime=acqTime, charge_state=header[0]['charge'],
                        precursor_mz=header[0]['precmz'], start_time=header[0]['rt'],
                        survey_id=int(header[0]['survey_spec']), parent_ion=header[0]['setmass'],
                        peak_rt=header[0]['rtapex'], s2i=header[0]['s2i'], peak_fwhm=header[0]['fwhm'],
                        fragenergy=header[0]['fragenergy'], precursor_noise=header[0]['thresh'],
                        sumreporterarea=header[0]['sumreparea'], sumreporterint=header[0]['sumrepint'],
                        peak_inten_raw=header[0]['inten'], peak_inten_smoothed=header[0]['smooth'],
                        peak_inten_integrated=header[0]['integ'], peak_area=header[0]['area'],
                        precursor_area=header[0]['c12'])

        specDict['sample_id'] = sample_id

        return specDict

    def makePeptideDictList(self, spectrum_id, pep, query, sets, hdf, pep2unique):

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
                       is_reverse_hit=0, missed_cleavage_sites=pep['missed_cleavage_sites'],
                       delta_seq=query['delta_seq'], ppm_error=pep['da_delta'] / pep['mw'] * 1e6,
                       delta_mod=query['delta_mod'], in_protein_inference=0,
                       fdr_at_score=score2fdr[round(pep['score'])][3])
        mods = hdf.modsHandler.generateModifications(pep)
        quantFound = 0
        # here we check that the modification relating to the quantification label is in the
        # mascot-assigned modifications
        if quantMethID:
            if quantMethData['mascot_name'] in mods['positional_modstring']:
                quantFound = 1
                pepDict['is_quantified'] = 1

        pepDict.update(mods)
        pepDictList = []
        for acc in sets:
            pepSet = self.setsManager.peptideSets[acc]
            seq_start = pepSet.peptideData[pep['peptide']]['seq_start']
            seq_end = pepSet.peptideData[pep['peptide']]['seq_end']

            if pepSet.is_reverse_hit:
                pepDict['is_reverse_hit'] = 1
            if pepDict['peptide'] in pepSet.validpeptides:
                pepDict['in_protein_inference'] = 1
            pepDict['seq_start'] = seq_start
            pepDict['seq_end'] = seq_end
            pepSet.adddataFromSet(pepDict)
            pepDictList.append(pepDict.copy())

        return pepDictList, quantFound

    def makeQuantDictLists(self, spectrum_id, specID, pepData, headerref, quanArray, hdf):
        score2fdr = self.importData.score2fdr
        quant = self.filterArrayEqual(quanArray, 'spec_id', specID)
        # all data in this list should be based on the rank 1 peptide sequence. If there are more than one entry it's
        # either due to more than one proteinset membership or peptide of different rank but the same score and
        # sequence composition.
        pepData.sort(key=lambda x: x['rank'])
        if len(pepData) > 1:
            is_unique = 0
        else:
            # take lowest unique value as an approximation for true uniqueness.
            is_unique = min(p['is_unique'] for p in pepData)
        firstpep = pepData[0]
        score = round(firstpep['score'])
        fdr_at_score = score2fdr[score][3]
        delta_seq = firstpep['delta_seq']
        peptide_length = len(firstpep['peptide'])
        quantList = []

        protein_group_nos = set([x['protein_group_no'] for x in pepData])
        deltas = []
        if len(quant) > 0:
            # accumulate numbers
            s2i = headerref['s2i']
            p2t = headerref['c12'] / headerref['thresh']
            hdf.quantified_spectra += 1
            if len(protein_group_nos) > 1:
                logger.log.debug('more than one proteinset (%s) !! %s' % (str(protein_group_nos), specID))
            # add to the hdf5results list
            for q in quant:
                for protein_group_no in protein_group_nos:
                    specQuantDict = dict(protein_group_no=protein_group_no, spectrum_id=spectrum_id,
                                         delta_seq=delta_seq, is_unique=is_unique,
                                         isotopelabel_id=int(q['isolabel_id']), score=score,
                                         quant_raw=q['area'], quant_isocorrected=q['corrected'],
                                         quant_allcorrected=0, p2t=p2t,
                                         peptide_length=peptide_length, s2i=s2i,
                                         in_quantification_of_protein=0, fdr_at_score=fdr_at_score)

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
        logger.setProtInferanceLogs()

        if resultfile.exists():
            resultfile.unlink()
        expID = None

        quantMethID = cfg.parameters['runtime']['quantmethod_id']
        if quantMethID:
            quantMethData = QuantMethods().getMethodByID(quantMethID)
            quantSource = quantMethData['source']
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
            searchData = {}
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
