#!/usr/bin/env python

"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This is the main runner to perform internalization of Mascot results to .hdf5
files, protein inference and then protein quantification

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here:
https://github.com/protcode/isob/
"""

# python libraries
import sys
import os
import time
from pathlib import Path

# CommonUtils libraries

from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.LoggingManager import Logger
import CommonUtils.QuantMethodHandler as qmh
from CommonUtils.hdf5Base import hdf5Base
from CommonUtils.hdf5Mascot import HDF5Mascot
from CommonUtils.hdf5Results import HDF5Results
import CommonUtils.ExceptionHandler as ExHa


logFormat = ' --logging.logdir %s --logging.loglevel %s --logging.screenlevel %s'


def checkFilePresence(dataDir):
    logger.log.info('Checking presence of required files ...')

    look_for = "FILE"

    hdf5s = []
    foundFiles = {}
    missingFiles = {}
    fileTypes = {'.hdf5': hdf5s}
    for file in dataDir.iterdir():
        try:
            fileTypes[file.suffix].append(file.name)
        except KeyError:
            # file type not required
            pass

    for idx, datFile in enumerate(dataDir.glob('*.dat')):
        logger.log.log(logger.PROCESS, '\t - %s :', datFile.name)
        datName = datFile.name
        with open(str(datFile), "r") as file_to_read:
            for line in file_to_read:
                items = line.strip().replace('\\', '/').split('=')
                if look_for in items:
                    if items[1]:
                        fileName = Path(items[1]).with_suffix('.hdf5').name
                        file2find = dataDir.joinpath(fileName)
                        if file2find.exists():
                            logger.log.log(logger.PROCESS, '\t\thdf5 file found: %s' % file2find.name)
                            imports = getImports(file2find)
                            impType = 0
                            if imports:
                                logger.log.log(logger.PROCESS, '\t\thdf5 has Mascot imports: %s' % ','.join(imports))
                                if len(imports) == 1 and datName in imports:
                                    impType = 1
                                else:
                                    impType = 2
                            foundFiles[datName] = {'.hdf5': fileName, 'imports': impType}
                        else:
                            try:
                                missingFiles[datName].append(fileName)
                            except KeyError:
                                missingFiles[datName] = [fileName]
                    break
    return foundFiles, missingFiles


def runMascotParser(datFile, hdf5file, dataDir):

    logger.log.log(logger.PROCESS, 'Starting mascotParser for %s' % datFile)

    cwd = os.getcwd()
    os.chdir('mascotParser')
    general = cfg.parameters['general']
    cmdline = 'python mascotParser.py --datadir %s --datfile %s' + logFormat
    cmdline += ' --hdf5file %s --general.minhooklength %s --general.delta %s'
    cmdline = cmdline % (dataDir, datFile, cfg.parameters['logging']['logdir'], cfg.parameters['logging']['loglevel'],
                         cfg.parameters['logging']['screenlevel'], hdf5file, general['minhooklength'], general['delta'])

    logger.log.debug('cmd: %s' % cmdline)
    os.system(cmdline)
    os.chdir(cwd)

    hdf5path = Path(hdf5file).stem + '.error'
    ret = fetchError(dataDir.joinpath(hdf5path))
    ret['cmdline'] = cmdline

    return ret


def runProteinInference(fileList):

    fileListString = ','.join(fileList)
    timeString = time.strftime('%Y%m%d_%H%M', time.localtime(time.time()))
    suffix = 'results_%s.hdf5' % timeString
    if len(fileList) > 1:
        infoStr = 'Starting merged protein inference for %i files'
        resultFile = '%s_merged_%s' % (dataDir.name[:25], suffix)
        logger.log.log(logger.PROCESS, infoStr % len(fileList))
    else:
        infoStr = 'Starting protein inference for %s'
        resultFile = '%s_%s_%s' % (dataDir.name[:25],  Path(fileList[0]).stem, suffix)
        logger.log.log(logger.PROCESS, infoStr % fileListString)

    cwd = os.getcwd()
    os.chdir('proteinInference')

    general = cfg.parameters['general']
    cmdline = 'python proteinInference.py --datadir %s --filelist %s --resultfilename %s' + logFormat
    cmdline += ' --general.fdrthreshold %s '
    cmdline = cmdline % (str(dataDir), fileListString, resultFile, cfg.parameters['logging']['logdir'],
                         cfg.parameters['logging']['loglevel'], cfg.parameters['logging']['screenlevel'],
                         general['fdrthreshold'])

    if quantMeth:
        cmdline += ' --quantmethod_id %s' % quantMeth['meth_id']
    else:
        cmdline += ' --quantmethod_id none'

    logger.log.debug('cmd: %s' % cmdline)
    os.system(cmdline)

    resultPath = dataDir.joinpath(resultFile)
    ret = fetchError(dataDir.joinpath(resultPath.stem + '.error'))
    ret['resultfile'] = resultPath
    ret['suffix'] = suffix
    ret['cmdline'] = cmdline
    if ret['code'] == 0:

        # add in the parameters from pymssafe if successfully ran protein inference
        os.chdir(cwd)
        preCfg = ConfigManager(Path('./preMascot.cfg'))
        hdf = HDF5Results(resultPath)
        hdf.appendOpen()
        hdf.addConfigParameters(preCfg.parameters, 'preMascot', '1 pyMSsafe')
        hdf.close()

        ret = fetchError(dataDir.joinpath(resultPath.stem + '.error'))
        ret['resultfile'] = resultPath
        ret['suffix'] = suffix
        ret['cmdline'] = cmdline

    os.chdir(cwd)
    return ret


def runProteinQuantification(resultsFile, suffix):
    cwd = os.getcwd()
    os.chdir('proteinQuantification')

    if cfg.parameters['general']['run_corrects2iquant']:
        logger.log.log(logger.PROCESS, 'running correctS2Iquant for %s' % resultsFile.name)
        cmdline = 'python correctS2Iquant.py --datadir %s --filename %s' + logFormat
        cmdline = cmdline % (str(dataDir), resultsFile.name, cfg.parameters['logging']['logdir'],
                             cfg.parameters['logging']['loglevel'], cfg.parameters['logging']['screenlevel'])
        logger.log.debug('cmd: %s' % cmdline)
        os.system(cmdline)

        # test for errors in correctS2Iquant
        resultPath = dataDir.joinpath(resultsFile)
        ret = fetchError(dataDir.joinpath(resultPath.stem + '.error'))
        ret['resultfile'] = resultPath
        ret['suffix'] = suffix
        ret['cmdline'] = cmdline
    else:
        ret = dict(code=0)

    if ret['code'] == 0:
        logger.log.log(logger.PROCESS, 'running protein quant for %s' % resultsFile.name)
        cmdline = 'python quantifyProteins.py --datadir %s --filename %s' + logFormat
        cmdline += ' --reference %i --general.fdrthreshold %s'
        cmdline = cmdline % (str(dataDir), resultsFile.name, cfg.parameters['logging']['logdir'],
                             cfg.parameters['logging']['loglevel'], cfg.parameters['logging']['screenlevel'],
                             refID, cfg.parameters['general']['fdrthreshold'])

        logger.log.debug('cmd: %s' % cmdline)
        os.system(cmdline)
        os.chdir(cwd)

        resultPath = dataDir.joinpath(resultsFile)
        ret = fetchError(dataDir.joinpath(resultPath.stem + '.error'))
        ret['resultfile'] = resultPath
        ret['suffix'] = suffix
        ret['cmdline'] = cmdline

    os.chdir(cwd)
    return ret


def runOutput(suffix):
    logger.log.log(logger.PROCESS, 'creating output for files like %s' % suffix)

    cmdline = 'python outputResults.py --datadir %s --filefilter *%s' + logFormat
    cmdline = cmdline % (str(dataDir), suffix, cfg.parameters['logging']['logdir'],
                         cfg.parameters['logging']['loglevel'], cfg.parameters['logging']['screenlevel'])
    os.system(cmdline)

    ret = fetchError(dataDir.joinpath(Path(suffix).stem + '.error'))
    ret['suffix'] = suffix
    ret['cmdline'] = cmdline

    return ret


def fetchError(errorFile):
    if errorFile.exists():
        # check to see if named file error exists: report
        fin = open(str(errorFile), 'r')
        error = eval(fin.readline())
        fin.close()
        return error
    else:
        genErrorFile = errorFile.parent.joinpath('errors.error')
        if genErrorFile.exists():
            # check to see if general error file exists: rename and report
            fin = open(str(genErrorFile), 'r')
            error = eval(fin.readline())
            fin.close()
            return error
        else:
            return dict(code=0, error=None)


def getQuantMethod(hdf5File):
    # reads the quantification method from the isotopes table.
    hdfMascot = HDF5Mascot(hdfFilePath=dataDir.joinpath(hdf5File))
    hdfMascot.readOpen()
    methodID = hdfMascot.getH5DFQuantMeth()
    hdfMascot.close()
    quantMethods = qmh.QuantMethods()
    methodData = quantMethods.getMethodByID(methodID)
    return methodData


def assignReferneceLabel():
    # assigns the referenceID for quantification, converting first and last to the lowest and highest
    # mass isotope respectively
    referenceID = None

    reference = cfg.parameters['quantification']['reference'].lower()
    labels = []
    qMasses = quantMeth['quantmasses']
    for id in qMasses:
        labels.append((id, qMasses[id][0]['mass'], qMasses[id][0]['name'].lower()))

    labels.sort(key=lambda x: x[1])

    if reference == 'first':
        referenceID = labels[0][0]
    elif reference == 'last':
        referenceID = labels[-1][0]
    else:
        for lab in labels:
            if reference == lab[2]:
                referenceID = lab[0]
            elif reference.isdigit() and int(reference) == lab[0]:
                referenceID = lab[0]

            if referenceID:
                break
    return referenceID


def getImports(hdf5File):
    # retrieves any previous Mascot imports
    hdfMascot = HDF5Mascot(hdfFilePath=hdf5File)
    hdfMascot.readOpen()
    try:
        imports = hdfMascot.hdf.readTable('/imports')
    except ExHa.HDF5consistancyError:
        hdfMascot.close()
        return []

    hdfMascot.close()
    if len(imports) == 0:
        return []
    else:
        return [x['datfile'] for x in imports]


def addIsotopeLabelTable(msHdf5, ResultHdf5):
    # transfers the isotopes table from the MS.hdf5 file to the results.hdf5 file.
    msHdfObj = hdf5Base(dataDir.joinpath(msHdf5))
    msHdfObj.readOpen()
    isotopeTable = msHdfObj.readTable('/rawdata/isotopes')
    msHdfObj.close()

    resultsHdfObj = hdf5Base(ResultHdf5)
    resultsHdfObj.appendOpen()
    resultsHdfObj.createTable('/', 'isotopes', 'IsotopeLabel', True)
    resultsHdfObj.appendRows('/isotopes', isotopeTable)
    resultsHdfObj.close()


def askQuestionAbortDeleteUse():
    # all files have mascot data so need to ask abort / delete / use existing
    userResponse = 0
    runMascot = True
    responses = ['1', '2', '3']
    while userResponse not in responses:
        print '\nMascot data has been detected in all .HDF5 files.  Do you want to?'
        print '  1) Abort the analysis'
        print '  2) Delete the existing Mascot data and use the new Mascot data for the analysis'
        print '  3) Continue the analysis with the current data'
        userResponse = raw_input('Please enter 1, 2 or 3: ')

        if userResponse in responses:
            if userResponse == '1':
                logger.log.warning('Aborting analysis')
                sys.exit()
            elif userResponse == '2':
                logger.log.log(logger.PROCESS, 'Deleting existing Mascot data')
                runMascot = True
            elif userResponse == '3':
                logger.log.log(logger.PROCESS, 'Using existing Mascot data')
                runMascot = False
        else:
            print '\nInvalid entry, please try again'
    return runMascot


def askQuestionAbortDelete():
    # only some files have mascot data so need to ask abort / delete
    userResponse = 0
    runMascot = True
    responses = ['1', '2']
    while userResponse not in responses:
        print '\nMascot data has been detected in the .HDF5 file.  Do you want to?'
        print '  1) Abort the analysis'
        print '  2) Delete the existing Mascot data and use the new Mascot data for the analysis'
        userResponse = raw_input('Please enter 1 or 2: ')

        if userResponse in responses:
            if userResponse == '1':
                logger.log.warning('Aborting analysis')
                sys.exit()
            elif userResponse == '2':
                logger.log.log(logger.PROCESS, 'Deleting existing Mascot data')
                runMascot = True
        else:
            print '\nInvalid entry, please try again'
    return runMascot

if __name__ == '__main__':

    configPath = './postMascot.cfg'
    cfg = ConfigManager(configPath)
    ret = cfg.evaluateCommandLineArgs(sys.argv)

    pID = os.getpid()

    dataDir = cfg.parameters['runtime']['datadir']
    logParam = cfg.parameters['logging']
    logPath = Path(dataDir.joinpath(logParam['logdir']))
    if logPath.exists():
        # clear existing logs for this part of the pipeline
        for log in logPath.glob('*.log'):
            if log.name in ['postMascot.log', 'mascotParser.log', 'proteinQuant.log',
                            'correctS2Iquant.log', 'output.log']:
                log.unlink()
    else:
        logPath.mkdir(parents=True)
    logFile = logPath.joinpath(logParam['logfile'])

    logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], pID)
    logger.setLog('postMascot')

    logger.log.log(logger.PROCESS, 'started postMascot application')
    logger.log.info('Removing previous errors from %s' % str(dataDir))
    # clear previous error files
    for errorFile in dataDir.glob('*.error'):
        errorFile.unlink()
    logger.log.info('Analysing: %s' % str(dataDir))

    foundFiles, missingFiles = checkFilePresence(dataDir)

    startProcessing = False
    if not foundFiles:
        logger.log.warning('No *.dat files found for processing. Please make sure these are in the same directory '
                           'as the .hdf5 files from the preMascot workflow.\nAborting.')

    elif missingFiles:
        if foundFiles:
            print 'Not all required files could be found:'
            for mf in missingFiles:
                print '\t%s : %s not found' % (mf, ' , '.join(missingFiles[mf]))

            userResponse = None
            while userResponse not in ['Y', 'N']:
                userResponse = raw_input('Start processing only for the present files? [Y/N]: ').upper()

            if userResponse == 'Y':
                startProcessing = True

        else:
            logger.log.warning('No files found for processing.\nAborting.')

    else:
        logger.log.info('All required files found.')
        startProcessing = True
    datList = sorted(foundFiles.keys())
    runMascot = False
    quantMeth = None
    if startProcessing:
        filesWimports = [foundFiles[x]['imports'] for x in foundFiles if foundFiles[x]['imports']]
        if len(filesWimports) == len(foundFiles):
            # all files have imports
            if max(filesWimports) > 1:
                # some files have multiple imports so cant use existing data
                runMascot = askQuestionAbortDelete()
            else:
                runMascot = askQuestionAbortDeleteUse()
        elif len(filesWimports) == 0:
            # no files have imports
            runMascot = True
        elif len(filesWimports) != len(foundFiles):
            # different numbers of files have imports
            runMascot = askQuestionAbortDelete()

        files2merge = []

        for datFile in datList:
            hdfFile = foundFiles[datFile]['.hdf5']
            quantMeth = getQuantMethod(hdfFile)
            if quantMeth:
                refID = assignReferneceLabel()

            if runMascot:

                ret = runMascotParser(datFile, hdfFile, dataDir)
                if ret['code'] == 0:
                    logger.log.info('Finished MascotParser successfully for %s' % hdfFile)
                else:
                    logger.log.warning('MascotParser failed for %s: %s' % (hdfFile, ret['error']))
                    logger.log.warning('cmd: %s' % ret['cmdline'])
                    continue

            if cfg.parameters['runtime']['mergeresults']:
                files2merge.append(hdfFile)
            else:
                ret = runProteinInference([hdfFile])
                if ret['code'] == 0:
                    logger.log.info('Finished proteinInference successfully for %s' % hdfFile)
                else:
                    logger.log.warning('proteinInference failed for %s: %s' % (hdfFile, ret['error']))
                    logger.log.warning('cmd: %s' % ret['cmdline'])
                    continue

                if quantMeth:
                    addIsotopeLabelTable(hdfFile, ret['resultfile'])
                    if ret['resultfile'].exists():
                        ret = runProteinQuantification(ret['resultfile'], ret['suffix'])

                        if ret['code'] == 0:
                            logger.log.info('Finished proteinQuantification successfully for %s' % hdfFile)
                        else:
                            logger.log.warning('proteinQuantification failed for %s: %s' % (hdfFile, ret['error']))
                            logger.log.warning('cmd: %s' % ret['cmdline'])
                            continue
                    else:
                        logger.log.warning('Could not find results file for %s' % ret['resultfile'])

                runOutput(ret['suffix'])

        if files2merge:
            ret = runProteinInference(files2merge)
            if ret['code'] == 0:
                logger.log.info('Finished proteinInference successfully for %s' % ret['resultfile'])
            else:
                logger.log.warning('proteinInference failed for %s: %s' % (ret['resultfile'], ret['error']))
                logger.log.warning('cmd: %s' % ret['cmdline'])
                sys.exit()

            if quantMeth:
                addIsotopeLabelTable(files2merge[0], ret['resultfile'])
                if ret['resultfile'].exists():
                    ret = runProteinQuantification(ret['resultfile'], ret['suffix'])
                    if ret['code'] == 0:
                        logger.log.info('Finished proteinQuantification successfully for %s' % ret['resultfile'])
                    else:
                        logger.log.warning('proteinQuantification failed for %s: %s' %
                                           (ret['resultfile'], ret['error']))
                        logger.log.warning('cmd: %s' % ret['cmdline'])
                        sys.exit()
                else:
                    logger.log.warning('Could not find results file for %s' % ret['resultfile'])

            ret = runOutput(ret['suffix'])

        if ret['code'] == 0:
            logger.log.log(logger.PROCESS, 'Finished processing %i files' % len(foundFiles))
        else:
            logger.log.warning('Errors detected in processing.')
