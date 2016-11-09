__author__ = 'holfra'

from hdf5Base import hdf5Base
import numpy as np
from hdf5Tables import Config


class HDF5Config:
    def __init__(self, hdfFilePath=None, hdfBaseObj=None):

        if not hdfFilePath and not hdfBaseObj:
            raise (ValueError, "Neither hdfFilePath nor hdfBaseObj provided.")

        elif hdfFilePath and hdfBaseObj:
            if hdfBaseObj.filePath == hdfFilePath:
                self.hdf = hdfBaseObj
                self.hdfFilePath = hdfFilePath

            else:
                raise (ValueError, "The provided hdfFilePath and hdfBaseObj do not reference the same file.")

        elif hdfBaseObj and not hdfFilePath:
            self.hdf = hdfBaseObj
            self.hdfFilePath = hdfBaseObj.filePath

        elif hdfFilePath and not hdfBaseObj:
            self.hdfFilePath = hdfFilePath

            self.hdf = hdf5Base(hdfFilePath)

    def appendOpen(self):
        self.hdf.appendOpen()

    def readOpen(self):
        self.hdf.readOpen()

    def close(self):
        self.hdf.close()

    def writeMS1QuantConfig(self, cfg):
        # first clear out any old MS1 quant data
        self.hdf.removeRows('/rawdata/config', {'set': 'SILAC'})
        self.hdf.removeRows('/rawdata/config', {'set': 'MS1Quant'})

        # now add the new parameters
        rowsToAdd = []
        for param in cfg.__dict__:
            if param == 'connectstring':
                continue
            rowsToAdd.append({'set': 'MS1Quant',
                              'parameter': str(param),
                              'value': str(getattr(cfg, param))})

        self.hdf.appendRows('/rawdata/config', rowsToAdd)
