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

This file controls the preMascot workflow.  It identifies the files to
process and calls the applications required to process the files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/cellzome/isobarquant
"""

# python libraries
import sys
import os

from pathlib import Path


# CommonUtils libraries
sys.path.insert(0, '..')
from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.LoggingManager import Logger

# application libraries


def runPyMSsafe():
    cmdline = 'python "%s/pyMSsafe.py" %s' % (pyMSsafeDir, formatArgumentString(rawFile, 0))
    logger.log.debug(cmdline)
    os.chdir(pyMSsafeAppDir)
    os.system(cmdline)

    ret = fetchError(rawFile.parent.joinpath(rawFile.stem + '.error'))
    ret['cmdline'] = cmdline
    return ret


def runMGFcreation():
    cmdline = 'python "%s/mgf.py" %s' % (pyMSsafeDir, formatArgumentString(rawFile, True))
    logger.log.debug(cmdline)
    os.chdir(pyMSsafeAppDir)
    os.system(cmdline)

    ret = fetchError(rawFile.parent.joinpath(rawFile.stem + '.error'))
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


def formatArgumentString(filePath, isMGF):
    argStr = ''

    # generate file and directory information
    direc = str(filePath.parent).replace('\\', '/')
    if direc.find(' ') != -1:
        direc = '"%s"' % direc
    file = filePath.name
    hdf5 = file[:-3] + 'hdf5'
    if file.find(' ') != -1:
        file = '"%s"' % file
        hdf5 = '"%s"' % hdf5

    # generate logging information
    logFormat = ' --logging.logdir %(logdir)s --logging.loglevel %(loglevel)s --logging.screenlevel %(screenlevel)s'
    logging = cfg.parameters['logging']
    logStr = logFormat % logging

    general = cfg.parameters['general']
    if isMGF:
        # provide commandline for MGF scripts
        argStr = '--datadir %s --filefilter %s* --msmsfilters.tolppm %s --msmsfilters.tolmda %s' + logStr
        argStr = argStr % (direc, hdf5, general['tolppm'], general['tolmda'])
    else:
        # provide command line for pyMSsafe scripts
        argStr = '--datadir %s --filefilter %s* --pid %i --quant %s --general.tolppm %s --general.tolmda %s '
        argStr += '--xic.beforepeak %s --xic.afterpeak %s' + logStr
        argStr = argStr % (direc, file, pID, cfg.parameters['runtime']['quant'], general['tolppm'], general['tolmda'],
                           general['beforepeak'], general['afterpeak'])
    return argStr


if __name__ == '__main__':
    sysPlatform = sys.platform
    if sysPlatform != 'win32':
        raise OSError('Invalid operating system (%s) for running analysis, "win32" needed.' % sysPlatform)

    configPath = './preMascot.cfg'
    cfg = ConfigManager(configPath)
    ret = cfg.evaluateCommandLineArgs(sys.argv)

    pID = os.getpid()
    # instrument = cfg.parameters['runtime']['instrument']

    dataDir = cfg.parameters['runtime']['datadir']
    pyMSsafeAppDir = str(Path('./pyMSsafe').resolve())

    logParam = cfg.parameters['logging']
    logPath = Path(dataDir.joinpath(logParam['logdir']))
    if logPath.exists():
        # clear existing logs for this part of the pipeline
        for log in logPath.glob('*.log'):
            if log.name in ['preMascot.log', 'pyMSsafe.log', 'mgfCreation.log']:
                log.unlink()
    else:
        logPath.mkdir(parents=True)
    logFile = logPath.joinpath(logParam['logfile'])

    logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'])  # , pID)
    logger.setLog('preMascot')

    logger.log.log(logger.PROCESS, 'started preMascot application')
    logger.log.info('Removing previous errors from %s' % str(dataDir))
    # clear previous error files
    for errorFile in dataDir.glob('*.error'):
        errorFile.unlink()
    logger.log.info('Analysing: %s' % str(dataDir))

    rawFileList = []
    for idx, raw in enumerate(dataDir.glob('*.raw')):
        rawFileList.append((raw, idx + 1))
    logger.log.info('Found %i files to process' % len(rawFileList))
    if not rawFileList:
        logger.log.critical('Found %i files to process' % len(rawFileList))
    pyMSsafeDir = Path('./pyMSsafe').absolute()
    ret = dict(code=1)
    cwd = os.getcwd()
    for rawFile, idx in rawFileList:
        print
        logger.log.info('---- NEW FILE ----')
        logger.log.log(logger.PROCESS, 'starting %i/%i: processing %s' % (idx, len(rawFileList), rawFile.name))

        # run pyMSsafe process
        ret = runPyMSsafe()

        if ret['code'] == 0:
            logger.log.info('Finished pyMSsafe succesfully for %s' % rawFile.name)

            # run mgf creation process
            ret = runMGFcreation()
            if ret['code'] == 0:
                logger.log.info('Finished MGFcreation succesfully for %s' % rawFile.name)
            else:
                logger.log.warning('MGFcreation failed for %s: %s' % (rawFile.name, ret['error']))
                logger.log.warning('cmd: %s' % ret['cmdline'])
        else:
            logger.log.warning('pyMSsafe failed for %s: %s' % (rawFile.name, ret['error']))
            logger.log.warning('cmd: %s' % ret['cmdline'])

        os.chdir(cwd)
        logger.log.info('finished %i/%i: processing %s' % (idx, len(rawFileList), rawFile.name))
        print

    if ret['code'] == 0:
        logger.log.log(logger.PROCESS, 'Finished processing %i files' % len(rawFileList))
    else:
        logger.log.warning('Errors detected in processing.')
