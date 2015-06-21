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

This file performs some simple mathematical and statistical operations on data.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

import numpy as np
import math


class Statistics():
    def __init__(self):
        pass

    @staticmethod
    def _toArray(inData):
        if type(inData) is np.ndarray:
            return inData
        else:
            return np.asarray(inData)

    def median(self, inData):

        data = self._toArray(inData)
        return np.median(data)

    def mean(self, inData):

        data = self._toArray(inData)
        return np.mean(data)

    def meanStd(self, inData):

        data = self._toArray(inData)
        if len(data) > 3:
            return np.mean(data), np.std(data)
        else:
            return np.mean(data), None


class Numerics:
    def __init__(self):
        pass

    @staticmethod
    def linearInterpolation(x, p0, p1):
        y = p0[1] + (p1[1] - p0[1]) / (p1[0] - p0[0]) * (x - p0[0])
        return y

    @staticmethod
    def roundDown(x, digits):
        fac = pow(10, digits)
        y = x * fac
        z = math.floor(y)
        return round(z / fac, digits)

    @staticmethod
    def roundUp(x, digits):
        fac = pow(10, digits)
        y = x * fac
        z = math.ceil(y)
        return round(z / fac, digits)
