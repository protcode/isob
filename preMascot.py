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

This file controls the preMascot workflow.  It identifies the files to
process and calls the applications required to process the files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""

# python libraries
import sys
import os
import multiprocessing

from pathlib import Path


# CommonUtils libraries
#$sys.path.insert(0, '..')
from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.LoggingManager import Logger
from CommonUtils.tools import MasterQueue

# application libraries


class RunPymssafe:
    def execute(self, rundata):
        rawFile, runDir, logger = rundata
        cmdline = 'python "pyMSsafe.py" %s' % formatArgumentString(rawFile, 0)
        logger.log.debug(cmdline)

        pyMSsafeAppDir = runDir / Path('pyMSsafe')

        os.chdir(str(pyMSsafeAppDir.resolve()))

        os.system(cmdline)
        ret = fetchError(rawFile.parent.joinpath(rawFile.stem + '.error'))
        ret['cmdline'] = cmdline
        # change back to runDir when done
        os.chdir(str(Path(runDir).resolve()))
        # we probably don't need this as the return codes are not kept
        ret['cmdline'] = cmdline
        return ret



class RunCreateMGF():
    def execute(self, rundata):
        rawFile, runDir, logger = rundata
        cmdline = 'python "mgf.py" %s' % formatArgumentString(rawFile, 1)
        logger.log.debug(cmdline)

        pyMSsafeAppDir = runDir / Path('pyMSsafe')
        os.chdir(str(pyMSsafeAppDir.resolve()))

        os.system(cmdline)

        ret = fetchError(rawFile.parent.joinpath(rawFile.stem + '.error'))
        # we probably don't need this as the retun codes are not kept
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
        argStr = '--datadir %s --filefilter %s* --msmsfilters.tolppm %s --msmsfilters.tolmda %s --general.printchargedata %s' \
                 + logStr
        argStr = argStr % (direc, hdf5, general['tolppm'], general['tolmda'], general['printchargedata'])
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
    RunPymssafeJobs = MasterQueue
    for rawFile, idx in rawFileList:
        print
        logger.log.info('---- NEW FILE ----')
        logger.log.log(logger.PROCESS, 'starting %i/%i: processing %s' % (idx, len(rawFileList), rawFile.name))
        RunPymssafeJobs.todoList.put_nowait((rawFile, os.path.dirname(os.path.abspath(__file__)), logger))
    available_threads = cfg.parameters['general']['thread_count']
    # if user supplies more threads than available CPUs, reset
    if available_threads > multiprocessing.cpu_count():
        available_threads = multiprocessing.cpu_count()
    my_threads = [RunPymssafeJobs() for i in range(available_threads)]


    for thread in my_threads:
        thread.add_job(RunPymssafe())
        thread.setDaemon(True)
        thread.start()

    for thread in my_threads:
        thread.join()
    #  will wait for all threads to finish before continuing:
    RunPymssafeJobs.todoList.join()

    RunCreateMGFJobs = MasterQueue
    for rawFile, idx in rawFileList:
        print
        logger.log.info('---- NEW FILE ----')
        logger.log.log(logger.PROCESS, 'mgf creation of %i/%i: processing %s' % (idx, len(rawFileList), rawFile.name))
        RunCreateMGFJobs.todoList.put_nowait((rawFile, os.path.dirname(os.path.abspath(__file__)), logger))

    available_threads = cfg.parameters['general']['thread_count']
    # if user supplies more threads than available CPUs, reset
    if available_threads > multiprocessing.cpu_count():
        available_threads = multiprocessing.cpu_count()
    my_mgf_threads = [RunCreateMGFJobs() for i in range(available_threads)]

    for thread in my_mgf_threads:
        thread.add_job(RunCreateMGF())
        thread.setDaemon(True)
        thread.start()

    for thread in my_mgf_threads:
        thread.join()
    RunCreateMGFJobs.todoList.join()

    logger.log.log(logger.PROCESS, 'Finished processing %i files' % len(rawFileList))

