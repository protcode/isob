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

This file mainly does processing of MS1 spectra to efficiently extract XIC
data

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""

import numpy as np
import math


class XIC:
    """
    @brief a class to handle data extraction from HDF5, data management and XIC processing
    """

    def __init__(self, hdf5object, cfg):
        """
        @brief creates the XIC object and links to hdf5 file
        @param hdf5object <hdf5Base object>: giving the data extraction methods needed to access XIC data
        @param cfg <Cfg object>: containing the controlling parameters for XIC processing
        @return:
        """
        self.hdf = hdf5object
        self.largeRTblock = 0
        self.BlockEdgeHigh = 0
        self.isotopeRTblock = 0

        self.XICparams = cfg.parameters['xic']

        rtn = self.hdf.getDataEqual('/rawdata/config', 'parameter', 'binsperda')
        if rtn:
            self.binsPerDa = int(rtn[0]['value'])
        else:
            self.binsPerDa = int(self.XICparams['binsperda'])

        surveyArray = hdf5object.getDataEqual('/rawdata/spectra', 'type', 'ms')
        self.surveyAll = surveyArray

        surveyDic = {}
        idxRT = {}
        surveyRT = {}
        for idx, s in enumerate(surveyArray):
            surveyDic[s['spec_id']] = idx
            idxRT[idx] = s['rt']
            surveyRT[s['spec_id']] = s['rt']

        self.maxIdx = idx
        self.surveyDic = surveyDic
        self.idxRT = idxRT
        self.surveyRT = surveyRT

    def setBlockEdge(self, RT, edge, width):
        """
        @brief finds the new high RT edge of the largeRTblock
        @param RT <float>: the new RT that doesn't fit to the old block
        @param edge <float>: current high RT edge of the current block
        @param width <float>: the block width
        @return:
        """
        newEdge = edge + width
        while RT > newEdge:
            newEdge += width

        return newEdge

    def loadLargeRTblock(self, RT):
        """
        @brief loads data from HDF5 file bounded by the RT limits
        @param RT <float>: upper RT limit of the block
        @return:
        """

        RToverlap = float(self.XICparams['rtoverlap'])
        blockWidth = float(self.XICparams['blockwidth'])
        RToverSmall = float(self.XICparams['rtoversmall'])
        blockWsmall = float(self.XICparams['blockwsmall'])

        tempEdgeHigh = self.setBlockEdge(RT, self.BlockEdgeHigh, blockWidth)

        lowRT = tempEdgeHigh - RToverlap - blockWidth
        highRT = tempEdgeHigh + RToverlap

        try:
            largeRTblock = self.hdf.getDataRange('/rawdata/xicbins', 'rt', lowRT, highRT)
        except MemoryError:
            # handle memory error
            tempEdgeHigh = self.setBlockEdge(highRT, self.BlockEdgeHigh, blockWsmall)
            lowRT = tempEdgeHigh - RToverSmall - blockWsmall
            highRT = tempEdgeHigh + RToverSmall
            largeRTblock = self.hdf.getDataRange('/rawdata/xicbins', 'rt', lowRT, highRT)

        self.BlockEdgeHigh = tempEdgeHigh
        self.largeRTblock = largeRTblock
        return

    def setXICblock(self, lowMZ, highMZ, lowRT, highRT):
        """
        @brief creates a new isotopBlock for XIC analysis
        @param lowMZ <float>: the low mass of XIC trace
        @param highMZ <float>: the high mass of XIC trace
        @param lowRT <float>: the low RT of the allowed peak range
        @param highRT <float>: the high RT of the allowed peak range
        @return:
        """

        binTol = int(self.XICparams['bintol_mda'])
        midRT = (highRT + lowRT) / 2

        if midRT > self.BlockEdgeHigh:
            # if the mid point of the new RT range is past the core of the current largeRTblock then load a new block
            self.loadLargeRTblock(midRT)

        # filter the largeRTblock to the desired mass range
        largeRTblock = self.largeRTblock
        lowBin = int(lowMZ * self.binsPerDa) - binTol
        highBin = int(highMZ * self.binsPerDa) + binTol
        xicFilter = (largeRTblock['bin'] >= lowBin) & (largeRTblock['bin'] <= highBin)
        self.isotopesBlock = largeRTblock[xicFilter]

        return

    def createXIC(self, mz, error):
        """
        @brief creates an XIC trace for a single mass with the given error
        @param mz <float>: mass for the center of the XIC
        @param error <float>: absolute error to define the range of the XIC
        @return:
        """

        isoBlock = self.isotopesBlock

        lowBin = int((mz - error) * self.binsPerDa)
        highBin = int((mz + error) * self.binsPerDa) + 1

        # self.xicData = isoBlock[(isoBlock['bin'] >= lowBin) & (isoBlock['bin'] <= highBin)]
        # return len(self.xicData)
        xicData = isoBlock[(isoBlock['bin'] >= lowBin) & (isoBlock['bin'] <= highBin)]

        prevSpec = -1
        prevIdx = -1
        removed = 0
        for idx in range(len(xicData)):
            if xicData[idx]['specid'] == prevSpec:
                # shared spectrum, select the closest to mz
                removed += 1
                prevIdx = self.selectXICionByMass(xicData, idx, prevIdx, mz)
                # prevIdx = self.selectXICionByIntensity(xicData, idx, prevIdx, mz)
                x = 1
            else:
                # spectrum different from last, nothing to do
                prevSpec = xicData[idx]['specid']
                prevIdx = idx

        xicData.sort(order='specid')
        self.xicData = xicData[removed:]
        return len(self.xicData)

    def selectXICionByMass(self, xicData, newIdx, oldIdx, mz):
        # shared spectrum, select the closest to mz
        deltaNew = abs(xicData[newIdx]['mz'] - mz)
        deltaOld = abs(xicData[oldIdx]['mz'] - mz)
        if deltaNew < deltaOld:
            # new is better
            xicData[oldIdx]['specid'] = 0
            oldIdx = newIdx
        elif deltaOld < deltaNew:
            # old is better
            xicData[newIdx]['specid'] = 0
        else:
            # offsets are the same, use the most intense
            if xicData[newIdx]['inten'] > xicData[oldIdx]['inten']:
                xicData[oldIdx]['specid'] = 0
                oldIdx = newIdx
            else:
                xicData[newIdx]['specid'] = 0

        return oldIdx

    def selectXICionByIntensity(self, xicData, newIdx, oldIdx, mz):
        if xicData[newIdx]['inten'] > xicData[oldIdx]['inten']:
            # new is better
            xicData[oldIdx]['specid'] = 0
            oldIdx = newIdx
        elif xicData[oldIdx]['inten'] > xicData[newIdx]['inten']:
            # old is better
            xicData[newIdx]['specid'] = 0
        else:
            # intensities the same select by mass
            if abs(xicData[newIdx]['mz'] - mz) < abs(xicData[oldIdx]['mz'] - mz):
                xicData[oldIdx]['specid'] = 0
                oldIdx = newIdx
            else:
                xicData[newIdx]['specid'] = 0

        return oldIdx

    def findMaxIntInTimeRange(self, lowRT, highRT, apex):
        xicData = self.xicData
        datalocation = (self.xicData['rt'] < highRT) & (self.xicData['rt'] > lowRT)
        apexlocation = xicData[datalocation]['rt'] == apex
        if apexlocation.any():
            apexInten = xicData[datalocation][apexlocation]['inten'][0]
            apexSpec = xicData[datalocation][apexlocation]['specid'][0]
        else:
            apexInten = 0
            apexSpec = 0

        surveydatalocation = (self.surveyAll['rt'] < highRT) & (self.surveyAll['rt'] > lowRT)
        specList = set(self.surveyAll[surveydatalocation]['spec_id'])
        newspeclist = specList - set(xicData[datalocation]['specid'])
        myapexSpec = 0
        for s in newspeclist:
            newspecapexloc = (self.surveyAll['spec_id'] == s) & (self.surveyAll['rt'] == apex)
            if newspecapexloc.any():
                myapexSpec = self.surveyAll[newspecapexloc]['spec_id'][0]

        if myapexSpec:
            apexSpec = myapexSpec
        if apexInten:
            return apexInten, apexSpec
        else:
            return 0, apexSpec

    def findPeaksInXIC(self, lowRT, highRT, offset):
        """
        @brief controls the peak detection from self.xicData, data is split into segments, processed and then peaks
        are found in the RT window
        @param lowRT <float>: low RT point of the RT window
        @param highRT <float>: high RT point of the RT window
        @param offset <integer>: offset from the monoisotopic ion
        @return:
        """
        xicData = self.xicData
        surveys = self.surveyDic
        lenData = len(xicData)
        maxGap = self.XICparams['clustgap']

        # initialise variables
        mzlow = 1000000
        mzhigh = 0
        first = 0
        last = 0
        end = 0
        peaks = []
        lastPosition = surveys[xicData[0]['specid']]

        for x in xrange(lenData):
            id = xicData[x]['specid']
            currentPosition = surveys[xicData[x]['specid']]

            if currentPosition - lastPosition > maxGap:
                end = x
            elif x == lenData - 1:
                end = lenData

            if end:
                if end - first > 1:
                    # more than one point to process
                    found = self.processSegment(xicData[first:end])
                    for pk in found:
                        # only append peaks in the RT range
                        if lowRT <= pk['rt'] <= highRT:
                            pk['offset'] = offset
                            peaks.append([pk])

                first = end
                end = 0

            lastPosition = currentPosition
            pass
        return peaks

    def findPeaksInXICNew(self, lowRT, highRT, offset):
        """
        @brief controls the peak detection from self.xicData, data is split into segments, processed and then peaks
        are found in the RT window
        @param lowRT <float>: low RT point of the RT window
        @param highRT <float>: high RT point of the RT window
        @param offset <integer>: offset from the monoisotopic ion
        @return:
        """
        xicData = self.xicData
        surveys = self.surveyDic
        lenData = len(xicData)
        maxGap = self.XICparams['clustgap']

        # initialise variables
        mzlow = 1000000
        mzhigh = 0
        first = 0
        last = 0
        end = 0
        peaks = []
        lastPosition = surveys[xicData[0]['specid']]

        for x in xrange(lenData):

            currentPosition = surveys[xicData[x]['specid']]

            if currentPosition - lastPosition > maxGap:
                end = x
            elif x == lenData - 1:
                end = lenData

            if end:
                if end - first > 1:
                    # more than one point to process
                    found = self.processSegment(xicData[first:end])
                    for pk in found:
                        # only append peaks in the RT range
                        if lowRT <= pk['rt'] <= highRT:
                            pk['offset'] = offset
                            peaks.append([pk])

                first = end
                end = 0

            lastPosition = currentPosition
            pass
        return peaks

    def processSegment(self, segment):
        """
        @brief performs peakwalking on the segment data, copying it into the smooth array first and performing smoothing
        @param segment <ndarray>: containing the raw ion data
        @return peaks <list>: containing the detected peaks in the segment
        """

        xicParams = self.XICparams
        calcRT50 = self.calcRT50
        processPeak = self.processPeak

        side = xicParams['av_window'] / 2
        mzlow = 100000
        mzhigh = 0

        # create new array to hold smoothed values
        smooth = np.ndarray(len(segment) + xicParams['av_window'] - 1,
                            dtype=[('spec', int), ('rt', float), ('mz', float),
                                   ('area', float), ('inten', float), ('smooth', float), ])

        # add blank data at the beginning and end
        row = side
        lowRT = self.idxRT[max(self.surveyDic[segment[0]['specid']] - 1, 0)]
        highRT = self.idxRT[min(self.surveyDic[segment[-1]['specid']] + 1, self.maxIdx)]


        for r in range(side):
            smooth.put(r, (0, lowRT, 0.0, 0.0, 0.0, 0.0))
            smooth.put(len(smooth) - 1 - r, (0, highRT, 0.0, 0.0, 0.0, 0.0))

        # copy seg data to the array
        for s in segment:
            smooth.put(row, (s['specid'], s['rt'], s['mz'], s['area'], s['inten'], 0.0))
            row += 1

        maxSmooth = 0
        minSmooth = 1
        localMaxPoint = 0
        localMinPoint = 0
        peaks = []
        rt50right = 0
        rt50left = 0
        scanRight = 0

        for point in xrange(side, len(smooth) - side):

            # calculate the new value
            # mean = np.average(smooth[point - side:point + side + 1])
            sum = 0.0
            for off in range(-side, side + 1):
                sum += smooth[point + off]['inten']

            mean = sum / xicParams['av_window']
            smooth[point]['smooth'] = mean

            # if below half height find the RT for the half height point
            if mean < maxSmooth * 0.5 and rt50right == 0:
                # find right half height
                rt50right = calcRT50(maxSmooth, smooth[point - 1], smooth[point])
                scanRight = point - 1

            if mean > maxSmooth:
                # increase the maximum and label as current best point
                maxSmooth = mean
                localMaxPoint = point
            elif mean < maxSmooth * xicParams['valley']:
                if maxSmooth / minSmooth > xicParams['threshold']:
                    peaks.append(
                        processPeak(smooth, localMaxPoint, localMinPoint, maxSmooth, point, rt50right, scanRight))

                minSmooth = mean
                maxSmooth = mean
                rt50right = 0
                localMaxPoint = point
                localMinPoint = point

            elif mean < minSmooth:
                if maxSmooth / minSmooth > xicParams['threshold']:
                    peaks.append(
                        processPeak(smooth, localMaxPoint, localMinPoint, maxSmooth, point, rt50right, scanRight))
                minSmooth = mean
                maxSmooth = mean
                rt50right = 0
                localMaxPoint = point
                localMinPoint = point
                pass

        if maxSmooth / minSmooth > xicParams['threshold']:
            peaks.append(processPeak(smooth, localMaxPoint, localMinPoint, maxSmooth, point + 1, rt50right, scanRight))

        return peaks

    def processPeak(self, smooth, apexPoint, minPoint, maxSmooth, point, rt50right, scanright):
        """
        @brief finds the edge points of the peak and calculates maxima and masses for the peak
        @param smooth <ndarray>: smoothed ion data
        @param apexPoint <integer>: index to the peak apex
        @param minPoint <integer>: index to the peak start point
        @param maxSmooth <float>: apex smoothed intensity
        @param point <integer>: index to the current location in the smooth array
        @param rt50right <float>: RT for the high RT 50% maximum point, if any
        @param scanright <integer>: index to the last location in the array greater than 50% maximum
        @return peak <dictionary>: containing the peak data
        """

        if not rt50right:
            rt50right = self.calcRT50(maxSmooth, smooth[point - 1], smooth[point])
            scanright = point - 1

        rt50left, scanleft = self.calcRT50left(smooth, maxSmooth, apexPoint, minPoint)
        # peakTop = smooth[scanleft:scanright + 1]
        mz = self.calcPeakMz(smooth, scanleft, scanright)

        peak = dict(specid=smooth[apexPoint]['spec'], rt=smooth[apexPoint]['rt'], smooth=smooth[apexPoint]['smooth'],
                    peakmz=mz['avg'], peaksd=mz['sd'], area=mz['maxarea'], inten=mz['maxint'], integ=mz['integ'],
                    rt50left=rt50left, rt50right=rt50right, scanleft=scanleft, scanright=scanright, apex=apexPoint,
                    c13=0)

        return peak

    def calcRT50(self, maxarea, pt1, pt2):
        """
        @brief interpolates the RT at 50% intensity point from the points either side of the crossover
        @param maxarea <float>: the peak maximum smooth value
        @param pt1 <dictionary>: point to the low RT side of the 50% crossover
        @param pt2 <dictionary>: point to the high RT side of the 50% crossover
        @return <float>: the RT of the 50% crossover point
        """
        return (maxarea * 0.5 - pt1['smooth']) * (pt2['rt'] - pt1['rt']) / (pt2['smooth'] - pt1['smooth']) + pt1['rt']

    def calcRT50left(self, smooth, maxSmooth, posApex, posMin):
        """
        @brief finds the low RT edge of the peak and calculates its RT
        @param smooth <ndarray>: smoothed ion data
        @param maxSmooth <float>: apex smoothed intensity
        @param posApex <integer>: index to the peak apex
        @param posMin <integer>: index to the peak start point
        @return rt <float>, pos <integer>: the RT of the 50% point and the index of the first point in the peak
        """

        found = 0
        halfHeight = maxSmooth / 2

        # find the first point below the halfHeight intensity
        for pos in xrange(posApex, posMin - 1, -1):
            if smooth[pos]['smooth'] < halfHeight:
                found = 1
                break

        if found:
            rt = self.calcRT50(maxSmooth, smooth[pos], smooth[pos + 1])
            pos += 1
        else:
            rt = smooth[posMin]['rt']
            pos = posMin
        return rt, pos

    def calcPeakMz(self, smooth, left, right):
        """
        @brief calculates m/z statistics of the peak apex
        @param peakTop <ndarray>: containing all points in the peak above 50% maximum smooth intensity
        @return mass <dictionary>: containing the calculated statistics
        """

        sumInt = 0.0
        sumProd = 0.0
        sumMZ = 0.0
        sumSq = 0.0
        integ = 0.0
        maxInten = 0.0
        maxArea = 0.0
        integ = 0.0
        n = right - left + 1

        for idx in range(left, right + 1):
            # for idx, ion in enumerate(peakTop):
            ion = smooth[idx]
            sumInt += ion['inten']
            sumProd += ion['inten'] * ion['mz']
            sumMZ += ion['mz']
            sumSq += ion['mz'] * ion['mz']

            if idx == left:
                integ += ion['inten'] * (smooth[idx + 1]['rt'] - ion['rt'])
            elif idx < right:
                integ += ion['inten'] * (smooth[idx + 1]['rt'] - smooth[idx - 1]['rt']) / 2
            else:
                integ += ion['inten'] * (ion['rt'] - smooth[idx - 1]['rt'])

            if ion['inten'] > maxInten:
                maxInten = ion['inten']

            if ion['area'] > maxArea:
                maxArea = ion['area']

        mass = dict(avg=sumProd / sumInt, mean=sumMZ / n, maxint=maxInten, maxarea=maxArea, integ=integ)

        if n == 1:
            mass['sd'] = 0.0
        else:
            numerator = (sumSq - sumMZ * sumMZ / n)
            if numerator < 1e-8:
                mass['sd'] = 0.0
            else:
                mass['sd'] = math.sqrt(numerator / n)

        return mass

    def calcPeakMz2(self, peakTop):
        """
        @brief calculates m/z statistics of the peak apex
        @param peakTop <ndarray>: containing all points in the peak above 50% maximum smooth intensity
        @return mass <dictionary>: containing the calculated statistics
        """

        sumInt = 0.0
        sumProd = 0.0
        sumMZ = 0.0
        sumSq = 0.0
        integ = 0.0
        maxInten = 0.0
        maxArea = 0.0
        n = len(peakTop)

        integ = peakTop[0]['inten'] * (peakTop[1]['rt'] - peakTop[0]['rt'])

        for idx, ion in enumerate(peakTop):
            sumInt += ion['inten']
            sumProd += ion['inten'] * ion['mz']
            sumMZ += ion['mz']
            sumSq = ion['mz'] * ion['mz']

            if idx < n - 1:
                integ += ion['inten'] * (peakTop[idx + 1]['rt'] - peakTop[idx - 1]['rt']) / 2
            else:
                integ += ion['inten'] * (ion['rt'] - peakTop[idx - 1]['rt'])

        mass = dict(avg=sumProd / sumInt, mean=sumMZ / n, maxint=maxInten, maxarea=maxArea, integ=integ)

        if n == 1:
            mass['sd'] = 0.0
        else:
            numerator = (sumSq - sumMZ * sumMZ / n)
            if numerator < 1e-8:
                mass['sd'] = 0.0
            else:
                mass['sd'] = math.sqrt(top / n)

        return mass

    def calcPeakMz3(self, peakTop):
        """
        @brief calculates m/z statistics of the peak apex
        @param peakTop <ndarray>: containing all points in the peak above 50% maximum smooth intensity
        @return mass <dictionary>: containing the calculated statistics
        """

        mass = dict(avg=np.average(peakTop['mz'], weights=peakTop['inten']),
                    mean=np.average(peakTop['mz']),
                    maxint=np.amax(peakTop['inten']),
                    maxarea=np.amax(peakTop['area']))

        if len(peakTop) == 1:
            integ = 0.0
            mass['sd'] = 0.0
        else:
            integ = peakTop[0]['inten'] * (peakTop[1]['rt'] - peakTop[0]['rt'])
            mass['sd'] = np.std(peakTop['mz'])

        points = len(peakTop)
        for pt in xrange(1, points):
            if pt < points - 1:
                integ += peakTop[pt]['inten'] * (peakTop[pt + 1]['rt'] - peakTop[pt - 1]['rt']) / 2
            else:
                integ += peakTop[pt]['inten'] * (peakTop[pt]['rt'] - peakTop[pt - 1]['rt'])

        mass['integ'] = integ

        return mass

    def linkOverlappingPeaks(self, basePeaks, newPeaks, isotope, rtOffset=0):

        for peakList in newPeaks:
            new = peakList[0]

            for baseList in basePeaks:
                base = baseList[0]
                prev = baseList[-1]

                if base['rt50left'] <= new['rt50right'] + rtOffset and base['rt50right'] >= new['rt50left'] - rtOffset:
                    # peaks overlap allowing for any rtOffset
                    # check to see if the last peak added was the same isotope as the previous addition
                    if prev['offset'] == new['offset']:
                        # exchange to the new peak if the mass difference to the theroretical is smaller
                        if abs(new['peakmz'] - isotope[1]) < abs(prev['peakmz'] - isotope[1]):
                            baseList[-1] = new.copy()
                    else:
                        baseList.append(new.copy())
                elif base['rt50left'] > new['rt50right'] + rtOffset:
                    # basePeak is past the range of the new peak so skip to the next peak
                    break
