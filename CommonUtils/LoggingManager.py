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

This file handles the logging interface.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

import logging
from pathlib import Path


class Logger:
    PROCESS = 45

    def __init__(self, logPath, logLevel='INFO', screenLevel='INFO', createNew=True):
        '''
        logfile needs to be a string and not a WindowsPath
        '''
        if createNew:
            fMode = 'w'
        else:
            fMode = 'a'

        try:
            logFile = str(logPath.absolute())
        except:
            logFile = str(Path(logPath).absolute())

        levels = {'DEBUG': logging.DEBUG,
                  'INFO': logging.INFO,
                  'WARNING': logging.WARNING,
                  'PROCESS': self.PROCESS}

        logging.basicConfig(level=levels[logLevel], format='%(asctime)s\t%(name)-12s\t%(levelname)-8s\t%(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S', filename=logFile, filemode=fMode)
        logging.addLevelName(self.PROCESS, 'PROCESS')

        # define a Handler for console output
        console = logging.StreamHandler()
        console.setLevel(levels[screenLevel])
        formatter = logging.Formatter('%(name)-14s: %(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    def setMascotParserLogs(self):
        self.setuplog = logging.getLogger('setup')
        self.datlog = logging.getLogger('datparse')
        self.qclog = logging.getLogger('dataQC')
        self.hdf5log = logging.getLogger('hdf5file')
        self.dblog = logging.getLogger('db')

    def setMSsafeLogs(self):
        self.log = logging.getLogger('pyMSsafe')
        self.dblog = logging.getLogger('db')

    def setProtInferanceLogs(self):
        self.log = logging.getLogger('protInference')
        self.hdf5log = logging.getLogger('hdf5log')

    def setLuxLogs(self):
        self.serverlog = logging.getLogger('Lux.server')
        self.queuelog = logging.getLogger('Lux.queue')
        self.logiclog = logging.getLogger('Lux.logic')

    def setLog(self, name):
        self.log = logging.getLogger(name)

    def setLogName(self, log, name):
        setattr(self, log, logging.getLogger(name))
