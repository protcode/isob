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

This file is the python interface to the czXcalibur.pyd library, which is the
C++ interface to the Xcalibur libraries.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

import czXcalibur as czX
import xcalparser as parser


class XRawFile(czX.XRawFile):
    """
    @brief class controlling access to the data in an Xcalibur RAW file
    @param czX.XRawFile <class>: this is the library class that XRawFile inherits properties from
    """

    def __init__(self, filename):
        """
        @Brief initialises the XRawFile object
        @param filename <string>: full path of the raw file
        """

        czX.XRawFile.__init__(self, filename)
        self.isopen = True
        self.__cache_numSpectra = None
        self.isExactive = 0
        self.skipScanEvents = []

    def close(self):
        """
        @Brief Close raw file
        """
        self.__cache_numSpectra = None
        self.isopen = False
        czX.XRawFile.close(self)

    def getNumSpectra(self):
        """
        @Brief Get number of spectra (scans)
        @Return 3-(tuple): (total number of spectra, number of first, number of last spectrum)
        """
        if not self.__cache_numSpectra:
            self.__cache_numSpectra = czX.XRawFile.getNumSpectra(self)
        return self.__cache_numSpectra

    def scanRange(self, firstScan=None, lastScan=None):
        """
        @Brief  helper method for looping over all existing (or subset) scans
        @Return <iterator> (xrange) of scan numbers
        """
        xfirst, xlast = self.getNumSpectra()[1:3]
        firstScan = max(firstScan, xfirst) if firstScan is not None else xfirst
        lastScan = min(lastScan, xlast) if lastScan is not None else xlast
        return xrange(firstScan, lastScan + 1)

    # And a property as syntactic sugar:
    scanNums = property(scanRange)

    def getAllFilters(self):
        """
        @Brief Get filter for all scans (spectra)
        @Return <list <string>> List of filters
        """
        return [self.getFilterForScanNum(scanNum) for scanNum in self.scanNums]

    # ##################### GavSwe additions

    def genSpecNum(self):
        """
        @brief generator function for the scan numbers
        """
        for scanNum in self.scanNums:
            yield scanNum

    def getSpectrumData(self, scanNum):
        """
        @brief fetches the spectrum ion data
        @param scanNum <integer>: spectrum number
        @return spec <list>: list of tuples containing (mz <float>, intensity <float>)
        """
        return self.getMassListFromScanNum(scanNum)

    def getSpectrumHeader(self, scanNum):
        """
        @brief fetches the spectrum ion data
        @param scanNum <integer>: spectrum number
        @return params <dictionary>: dictionary of all the parameters with their values.
        """
        params = dict(extra=parser.parseParameterList(self.getTrailerExtraForScanNum(scanNum)),
                      rt=self.getRTFromScanNum(scanNum),
                      filter=parser.parseFilter(self.getFilterForScanNum(scanNum))
                      )
        return params

    #    def getLCMethods(self, anal, load, dionex, pumps):
    def getLCMethods(self, methodsDic):
        """
        @brief retrieves the LC method parameters, designating analytical (fixed flow 'Qfixed') and loading LC methods
        @return lc_param <dictionary>: containing the analytical and loading LC method dictionaries of parameters
        """
        lc_param = dict(analytical={}, loading={}, lctype='Unknown')
        meths = self.getInstMethodNames()
        for j in range(len(meths)):
            if meths[j] in methodsDic['analytical']:
                tmp = parser.parseLCparameters(self.getInstMethod(j))
                tmp['MethodName'] = meths[j]
                lc_param['analytical'] = tmp.copy()
                lc_param['lctype'] = 'Eksigent'
            elif meths[j] in methodsDic['loading']:
                tmp = parser.parseLCparameters(self.getInstMethod(j))
                tmp['MethodName'] = meths[j]
                lc_param['loading'] = tmp.copy()
                lc_param['lctype'] = 'Eksigent'
            elif meths[j] in methodsDic['dionex']:
                lc_param = parser.parseDionexLCparameters(self.getInstMethod(j), methodsDic['pumps'])
                lc_param['lctype'] = 'Dionex'

        return lc_param

    def getMSMethods(self, meth):
        """
        @brief retrieves the MS method parameters
        @return ms_param <dictionary>: containing the MS method parameters
        """
        ms_param = {}
        order = []
        meths = self.getInstMethodNames()
        if 'Thermo Exactive' in meths:
            self.isExactive = 1
        for j in range(len(meths)):
            if meths[j] in meth:
                if self.isExactive:
                    # parsing for the Q Exactive
                    ms_param, order = parser.parseQExative(self.getInstMethod(j))
                else:
                    # parsing for the LTQ Orbitraps
                    ms_param, order = parser.parseLTQ(self.getInstMethod(j), self.skipScanEvents)

        if not self.isExactive:
            # find the mass range of the MS unit (Scan Event 1)
            if ms_param['Scan Event 1']['group'] == 'MS event':
                txt = ms_param['Scan Event 1']['Type']
                pos1 = txt.find('(')
                pos2 = txt.find('-')
                pos3 = txt.find(')')
                ms_param['Scan Event 1']['low mz'] = float(txt[pos1 + 1: pos2])
                ms_param['Scan Event 1']['high mz'] = float(txt[pos2 + 1: pos3])

            e1, e2, e3 = parser.parseEvents(ms_param, self.skipScanEvents)
            ms_param['units'] = e1
            ms_param['activation'] = e2
            ms_param['unit problems'] = e3

        return ms_param, order

    def getTuneParam(self):
        """
        @brief retrieves the MS tuning parameters
        @return tune_param <dictionary>: containing the MS tuning parameters
        """

        return parser.parseParameterList(self.getTuneData(0))
