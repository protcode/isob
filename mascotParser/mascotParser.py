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

This file is the main controller for the extraction of data from Mascot .dat
files and the storage of the data in MS.hdf5 files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""

# package imports
import sys
from pathlib import Path

# CommonUtils imports
sys.path.insert(0, '..')
from CommonUtils.tools import Stopwatch
from CommonUtils.hdf5Mascot import HDF5Mascot
from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.LoggingManager import Logger
import CommonUtils.ExceptionHandler as ExHa
from CommonUtils.QuantMethodHandler import QuantMethods

# project imports
from datfile import Datfile
from datparser import DatParser


def jobcontrol(cfg, logs):
    '''
    @brief takes filepath, executes parse ... import
    @param cfg <ConfigManager>: all the configuartion parameters including the sample data
    @param logs <loggingManager>: to control the logging of events
    '''
    try:
        msg = 'Setting files'
        datFileName = cfg.parameters['runtime']['datfile']
        hdf5FileName = cfg.parameters['runtime']['hdf5file']

        logs.setuplog.info('started job: %s' % datFileName)
        watch = Stopwatch()

        filedic = dict(dat=datFileName, datpath=dataDir.joinpath(datFileName),
                       hdf5=hdf5FileName, hdf5path=dataDir.joinpath(hdf5FileName))

        hdfMascot = HDF5Mascot(hdfFilePath=filedic['hdf5path'])
        hdfMascot.appendOpen()
        importGroup = hdfMascot.checkDatFilePresence(filedic['dat'])

        # test the quantification method already in the hdf5 file
        msg = 'Setting quantification method'
        quantMethID = hdfMascot.getH5DFQuantMeth()
        quantHandler = QuantMethods()
        quantMethod = quantHandler.getMethodByID(quantMethID)

    except Exception as genEx:
        ExHa.addContext(genEx, 'jobcontrol Error: %s' % msg)
        raise

    try:
        # control the deletion of existing data
        msg = 'Deleting existing HDF5 data'
        hdfMascot.deleteAllMascotImports(0)

        importGroup = filedic['dat'].replace('.', '_')

        hdfMascot.createTables(importGroup, 0)

        hdfMascot.writeConfig(cfg.convertConfig())

        datfile = Datfile(filedic, hdfMascot, cfg, logs, quantMethod)

        logs.setuplog.info('dat file: %s, hdf5 file: %s ' % (datfile.datfilename, filedic['hdf5']))

        msg = 'Parsing data'
        datparser = DatParser(datfile, cfg, logs)
        datparser.startParsing()
        watch.rec('parser')

        # post parsing processing
        logs.qclog.info('Effective Run-time analysis')
        datfile.doTimeAnalysis()

        msg = 'Find top protein hit'
        logs.qclog.info(msg)
        if datfile.seq2acc:
            datfile.findBestProtein()
    except Exception as genEx:
        ExHa.addContext(genEx, 'jobcontrol Error: %s' % msg)
        raise
    finally:
        hdfMascot.close()
        logs.setuplog.info('Closing HDF5')
        hdfMascot.close()
        watch.rec('processing')

        watch.stop()
        logs.setuplog.info('job took %s' % watch.format())
    return

# ########### MAIN ###############
if __name__ == '__main__':
    logs = 0
    cfg = ConfigManager('./mascotparser.cfg')
    ret = cfg.evaluateCommandLineArgs(sys.argv)

    try:
        cfg.scalePpmMda()
        dataDir = cfg.parameters['runtime']['datadir']

        logParam = cfg.parameters['logging']
        logPath = Path(dataDir.joinpath(logParam['logdir']))
        if not logPath.exists():
            logPath.mkdir(parents=True)
        logFile = logPath.joinpath(logParam['logfile'])

        logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], False)
        logger.setMascotParserLogs()

        jobcontrol(cfg, logger)

    except ExHa.UsageError as useEx:
        ExHa.reformatException(useEx)
        print useEx.context
    except Exception as genEx:
        ExHa.reformatException(genEx)
        errorFile = Path(cfg.parameters['runtime']['hdf5file']).stem + '.error'
        ExHa.exportError2File(genEx, cfg.parameters['runtime']['datadir'].joinpath(errorFile))
        if logs:
            logs.datlog.warning(ExHa.oneLineRepr(genEx))
        else:
            print ExHa.multiLineRepr(genEx)
