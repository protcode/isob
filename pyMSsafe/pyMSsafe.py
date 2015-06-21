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

This file controls the processing of Xcalibur .raw files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

__version__ = "$Revision: 1.0 $"

recpts = ['runtime', 'Load parameters', 'Process Spectra', 'XIC creation',
          'Writing MS/MS headers', 'Writing indexes']

# python modules
import xmlrpclib
import os
import sys
import socket
import traceback
from pathlib import Path

# cellzome CommonUtils
sys.path.insert(0, '..')
from CommonUtils.tools import *
from CommonUtils.ConfigManager import pyMSsafeConfigManager
from CommonUtils.LoggingManager import Logger
from CommonUtils.hdf5Base import hdf5Base
from CommonUtils.QuantMethodHandler import QuantMethods
import CommonUtils.ExceptionHandler as ExHa

# pyMSsafe modules
try:
    from xRawFile import XRawFile
except ImportError, ieEx:
    ExHa.reformatException(ieEx)
    ExHa.addContext(ieEx, 'Xcalibur not set up properly')
    configPath = './pymssafe.cfg'
    cfg = pyMSsafeConfigManager(configPath)
    ret = cfg.evaluateCommandLineArgs(sys.argv)
    dataDir = cfg.parameters['runtime']['datadir']
    ExHa.exportError2File(ieEx, dataDir.joinpath('errors.error'))

from datastore import Datamanager


class pymssafe:
    def __init__(self, config):
        """
        @brief initiaise pyMSsafe controller
        @param mode <string>: operational mode of the application
        """

        self.cfg = config

        dataDir = cfg.parameters['runtime']['datadir']
        logParam = cfg.parameters['logging']
        logPath = Path(dataDir.joinpath(logParam['logdir']))
        if not logPath.exists():
            logPath.mkdir(parents=True)
        logFile = logPath.joinpath(logParam['logfile'])
        self.logs = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], False)
        self.logs.setLog('pyMSsafe')

    def run(self, jobObj):
        """
        @brief run the analysis
        @param job <object>: containing all the job data
        """
        # needs to catch exceptions here so that LUX doesn't recieve an exception

        hdf = ''
        xraw = ''
        tempFile = ''
        config = self.cfg
        logs = self.logs

        try:
            hdf = ''
            maxspec = jobObj.args['maxspec']
            logs.log.info("Starting PyMSSafe runner")

            # assign objects to the file interfaces
            rawpath = jobObj.srcpth.absolute()
            namebase = rawpath.stem
            hdfpath = jobObj.dstpth.absolute().joinpath(namebase + '.hdf5')
            runname = namebase
            tempFile = jobObj.dstpth.absolute().joinpath(str(config.parameters['runtime']['pid']))
            if tempFile.exists():
                tempFile.unlink()

            logs.log.info('started job: %s' % rawpath.name)
            watch = Stopwatch()

            if not rawpath.exists():
                logs.log.info('could not find file: %s' % str(rawpath))
                logs.log.info('Stopped')
                return {'code': 1, 'error': 'could not find file: %s' % str(rawpath)}

            # get quant information if file is quan type
            quantMeth = self.extractQuant(jobObj.srcpth.name.upper())

            if quantMeth == -1:
                raise ExHa.QuantificationMethodError('Unable to find "%s" quantification method'
                                                     % self.cfg.parameters['runtime']['quant'].upper())
            elif quantMeth == -2:
                raise ExHa.QuantificationMethodError('Unable to find valid quantification method in file name (%s)'
                                                     % jobObj.srcpth.name.upper())

            xraw = XRawFile(str(rawpath))
            if config.parameters['general']['skipscanevents']:
                xraw.skipScanEvents = config.parameters['general']['skipscanevents']

            # opens hdf5 file for writing
            hdf = hdf5Base(hdfpath, True, True)
            hdf.appendOpen()
            self.createHDFtables(hdf)

            # save the config entries
            hdf.appendRows('/rawdata/config', config.convertConfig())

            # create datamanager object
            config.rawfile = jobObj.srcpth
            dataman = Datamanager(xraw, config, logs, hdf, quantMeth, str(tempFile))
            dataman.maxspec = maxspec
            dataman.addLUXdata(jobObj, config.parameters['runtime']['pid'])

            # run the analysis in the datamanager
            ok = dataman.processSpectra(quantMeth, watch, runname)
            if ok['code'] != 0:
                raise ExHa.SpectraProcessingError(ok['error'])

            # run the XIC generation
            logs.log.info('Processing XIC data for %d MS/MS events' % len(dataman.specs))
            ok = dataman.processXICs()
            if ok['code'] != 0:
                raise ExHa.XICprocessingError(ok['error'])
            watch.rec('Processing XIC')

            logs.log.info('Writing HDF5 indexes')
            hdf.indexTable('/rawdata/msmsheader', ['spec_id'])
            hdf.close()

            xraw = ''
            hdf = ''
            watch.stop()
            logs.log.info('job took: %s' % watch.oneLineFormat())
            tempFile.unlink()
        except ExHa.MSdataConsistancyError, msg:
            self.shutdown(hdf, xraw, tempFile)
            logs.log.warning('error: pyMSsafe Data error: %s' % msg)
            return {'code': 2, 'error': 'pyMSsafe Data error: %s' % msg}
        except ExHa.XICprocessingError, msg:
            self.shutdown(hdf, xraw, tempFile)
            logs.log.warning('error: pyMSsafe XIC error: %s' % msg)
            return {'code': 3, 'error': 'pyMSsafe XIC error: %s' % msg}
        except ExHa.SpectraProcessingError, msg:
            self.shutdown(hdf, xraw, tempFile)
            logs.log.warning('error: pyMSsafe Spectra error: %s' % msg)
            return {'code': 4, 'error': 'pyMSsafe Spectra error: %s' % msg}
        except ExHa.FragmentMethodError, msg:
            self.shutdown(hdf, xraw, tempFile)
            logs.log.warning('error: pyMSsafe Fragmentation method: %s' % msg)
            return {'code': 5, 'error': 'pyMSsafe Fragmentation method: %s' % msg}
        except ExHa.QuantificationMethodError, msg:
            self.shutdown(hdf, xraw, tempFile)
            logs.log.warning('error: pyMSsafe Quantification Batch Error: %s' % msg)
            return {'code': 7, 'error': 'Quantification Batch Error: %s' % msg}
        except ExHa.HDFwritingError, msg:
            self.shutdown(hdf, xraw, tempFile)
            logs.log.warning('error: pyMSsafe HDF5 error: %s' % msg)
            return {'code': 8, 'error': 'HDF5 error: %s' % msg}
        except:
            # other exceptions
            self.shutdown(hdf, xraw, tempFile)
            traceStr = repr(traceback.format_exception(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
            if logs != 'none':
                logs.log.warning('error: pyMSsafe General error: %s' % traceStr)
            return {'code': 1, 'error': 'pyMSsafe General error: %s' % traceStr}

        return {'code': 0}

    def extractQuant(self, fileName):
        qmh = QuantMethods()

        quantMethod = None
        cmdQuant = self.cfg.parameters['runtime']['quant'].upper()
        if cmdQuant == '0':
            # find quant in filename
            for name in qmh.methodsByName:
                if name in fileName:
                    quantMethod = qmh.methodsByName[name]
                    break
            if quantMethod is None:
                return None
        elif cmdQuant is None:
            return None
        else:
            # use method in file
            quantMethod = qmh.getMethodByName(cmdQuant)

        quantMethod['correction'] = qmh.extractCorrectionFactors(quantMethod)

        return quantMethod

    @staticmethod
    def shutdown(hdf, xraw, temp):
        try:
            if hdf:
                if hdf.isopen == 'r':
                    hdf.readclose()
                elif hdf.isopen == 'w':
                    hdf.writeclose()
            if xraw.isopen:
                xraw.close()

            if temp and temp.exists():
                temp.remove()

        except:
            pass

    @staticmethod
    def createHDFtables(hdf):
        """
        @brief calls the hdf module to create the appropriate tables and groups in the hdf5 file
        @param hdf <hdf5Base object>: hdf5Base module for manipulating hdf5 files
        @return:
        """

        hdf.createGroup('rawdata')

        tables = [('spectra', 'Spectra'), ('msmsheader', 'MSMSheader'), ('quan', 'Quan'), ('specparams', 'SpecParams'),
                  ('parameters', 'Parameters'), ('ions', 'Ions'), ('config', 'Config'), ('noise', 'Noise'),
                  ('isotopes', 'IsotopeLabel'), ('units', 'Units'), ('xicbins', 'XICbins')]
        for t in tables:
            hdf.createTable('rawdata', t[0], t[1])


# ########### MAIN ###############
if __name__ == '__main__':

    offline = 0
    configPath = './pymssafe.cfg'
    cfg = pyMSsafeConfigManager(configPath)
    ret = cfg.evaluateCommandLineArgs(sys.argv)
    cfg.scalePpmMda()
    pid = os.getpid()

    hostname = socket.gethostname().split('.')[0]
    client = "%s.%i" % (hostname, pid)

    mssafe = pymssafe(cfg)

    # default values for rawdir, rawfilter and inst
    rawdir = cfg.parameters['runtime']['datadir']

    mssafe.logs.log.info('Looking for: %s' % str(rawdir.joinpath(cfg.parameters['runtime']['filefilter'])))
    files = []
    for idx, raw in enumerate(rawdir.glob(cfg.parameters['runtime']['filefilter'])):
        files.append(Jobcontainer(raw, idx))

    if len(files) > 1:
        mssafe.logs.log.info('%d files to process' % len(files))

    i = 0
    for job in files:
        if job.srcpth.suffix.lower() != '.raw':
            mssafe.logs.log.log(mssafe.logs.PROCESS, 'Not RAW file, skipping: %s' % job.srcpth.name)
            continue

        i += 1
        if len(files) > 1:
            mssafe.logs.log.info('jobid = %s (%i/%i)' % (str(job.job_id), i, len(files)))
        # ret = dict(code=1, error='forced error')
        try:
            ret = mssafe.run(job)
        except:
            # handles exceptions in the exception handler of the controller
            trace = repr(traceback.format_exception(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
            ret = dict(code=99, error='Exception Handling Error: %s' % trace)

        if ret['code'] != 0:
            # error on the run
            errorFile = job.srcpth.parent.joinpath(job.srcpth.stem + '.error')
            fout = open(str(errorFile), 'w')
            fout.write(str(ret))
            fout.close()
