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

This file mainly provides processing methods for MS and MS/MS data.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/cellzome/isobarquant
"""

import math
import numpy
from pathlib import Path


class Processor():
    """
    @brief controlls the methods required to process profile spectrum data
    """

    def __init__(self, cfg):
        """
        @Brief Initiator of Processor class
        @Return None
        """
        self.smoothing = dict(mean=self.calcmean, median=self.calcmedian, savgol=self.savitzky_golay)
        self.pick = dict(centroid=self.centroid, centroid3=self.centroid3,
                         gaussian=self.gaussian, loggauss=self.loggaussian)
        self.cfg = cfg
        self.isMS2 = 0

        # self.isomodel = isoModel(cfg, 100, 0.01)

    def process(self, data, smooth='', pick='', filt='', isMS2=0):
        """
        @brief triggers the processing of the spectrum
        @param data <list>: contains the data for the spectrum processing
        @param smooth <string>: contains the type of smoothing required overwrites the config file method
        @param pick <string>: contains the type of peak picking required overwrites the config file method
        @param filt <string>: contains the type of ion filtering required overwrites the config file method
        @return filtered, mininten <
        """
        cfg = self.cfg
        self.isMS2 = isMS2

        if len(data) == 0:
            return [], 0

        # uses the config methods if none supplied
        if not smooth:
            smooth = cfg.parameters['smoothing']['smooth']
        if not pick:
            pick = cfg.parameters['picking']['pick']
        if not filt:
            filt = cfg.parameters['filtering']['filter']

        if smooth == 'savgol':
            processed = self.savitzky_golay(data)
        elif smooth:
            processed = self.smooth(data)
        else:
            processed = data[:]

        processed.insert(0, (processed[0][0] - 0.002, 0))
        processed.append((processed[-1][0] + 0.002, 0))

        # if pick:
        # perform the selected peak picking method
        pickmeth = self.pick[pick]

        # find the peak tops of all peaks incl noise!
        peaktops = self.findpeaktops(processed)

        # for each peak top find the limits of the peak
        for top in xrange(len(peaktops)):
            # integrate to find the area between the limits
            area = 0.0
            for point in xrange(peaktops[top]['left'], peaktops[top]['right'] + 1):
                area += processed[point][1] * (processed[point + 1][0] - processed[point - 1][0]) / 2
            peaktops[top]['area'] = area

            # find the accurate mass of the peak
            peaktops[top]['mz'] = pickmeth(peaktops[top], processed)
            peaktops[top]['logint'] = math.log10(peaktops[top]['topint'])
            if area > 0:
                peaktops[top]['logarea'] = math.log10(area)
            else:
                peaktops[top]['logarea'] = 0.0

        if filt == 'intensity':
            # perform filtering based on ion intensity
            filtered, mininten = self.filter4inten(peaktops)
        elif filt == 'charge':
            # perform filtering based on pairs of peaks matching to charge states
            filtered = self.filter4charge(peaktops)
            mininten = 0
        else:
            filtered = peaktops[:]
            for f in range(len(filtered)):
                filtered[f]['index'] = f
            mininten = 0

        return filtered, mininten

    def smooth(self, data):
        """
        @brief performs a mean smooth
        @params data <list>: spectrum data, list (mz <float>, inten <float>)
        @return smooth <list>: smoothed spectrum data, list (mz <float>, inten <float>)
        """
        resetblock = self.resetblock
        smoothParam = self.cfg.parameters['smoothing']
        smooth = self.smoothing[smoothParam['smooth']]

        dlen = len(data)
        window = smoothParam['smoothwindow']
        win2 = window / 2
        processed = []
        block = resetblock(data, 0, window)

        point = -1
        while point < dlen - 1:
            point += 1
            # add the next data point and remove the first
            block.pop(0)
            if point + win2 >= dlen:
                block.append((0.0, 0.0))
            else:
                block.append(data[point + win2])
            processed.append((data[point][0], smooth(block, window)))

            # if zerosep == 1 then two adjacent zero intensities indicate a break in the mass scale
            if smoothParam['zerosep']:
                if block[-1][1] == 0 and block[-2][1] == 0 and block[-1][0] > 0:
                    # last two data points have zero inten so end of peak
                    # process the end of the peak
                    while block[win2][1] > 0:
                        point += 1
                        block.pop(0)
                        block.append((0.0, 0.0))
                        processed.append((data[point][0], smooth(block, window)))

                    # now reset for the beginning of the next peak
                    block = resetblock(data, point + 1, window)

        return processed

    @staticmethod
    def resetblock(data, pos, window):
        """
        @brief determines the block of ions for smoothing
        @param data <list>: contains the MS spectrum data
        @param pos <integer>: positon for the block to start
        @params window <integer>: the number of data points in the block
        @return block <list>: the selected ions from the spectrum
        """
        block = [(0.0, 0.0)] * window

        # set starting positon
        for point in range(window / 2):
            block.pop(0)
            block.append(data[point + pos])
        return block

    @staticmethod
    def calcmean(data, window):
        """
        @brief calculates the average intensity of the block data
        @param data <list>: contains the MS spectrum data to be averaged
        @return avg <float>: the average of the values
        """
        avg = 0.0
        for point in data:
            avg += point[1]
        avg /= window
        return avg

    # noinspection PyUnusedLocal
    @staticmethod
    def calcmedian(data, window):
        """
        @brief calculates the median intensity of the block data
        @param data <list>: contains the MS spectrum data to find the mean of
        @return avg <float>: the median of the values
        """
        l = [x[1] for x in data]
        l.sort()
        num = len(l)
        mid = num / 2
        if num % 2 == 0:
            # even length use avge of center pts
            median = (l[mid] + l[mid - 1]) / 2
        else:
            # odd length use middle point
            median = l[mid]
        return median

    def savitzky_golay(self, data):
        """
        @brief applies a Savitzky-Golay filter
        @param data <list>: spectrum data as list of (mz <float>, inten <float>)
        @return smoothed data as a numpy array
        """
        try:
            kernel = self.cfg.parameters['smoothing']['smoothwindow']
            order = self.cfg.parameters['smoothing']['order']
        except ValueError:
            raise ValueError("kernel and order have to be of type int (floats will be converted).")
        if kernel % 2 != 1 or kernel < 1:
            raise ValueError("kernel size must be a positive odd number, was: %d" % kernel)
        if kernel < order + 2:
            raise ValueError("kernel is to small for the polynomals\nshould be > order + 2")

        # a second order polynomal has 3 coefficients
        order_range = range(order + 1)
        half_window = (kernel - 1) // 2
        b = numpy.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
        # since we don't want the derivative, else choose [1] or [2], respectively
        m = numpy.linalg.pinv(b).A[0]
        window_size = len(m)
        half_window = (window_size - 1) // 2

        # precompute the offset values for better performance
        offsets = range(-half_window, half_window + 1)
        offset_data = zip(offsets, m)

        smooth_data = []

        # temporary data, with padded zeros (since we want the same length after smoothing)
        pad = [(0.0, 0.0)] * half_window
        paddata = pad + data + pad

        # data = numpy.concatenate((numpy.zeros(half_window), data, numpy.zeros(half_window)))
        for i in range(half_window, len(paddata) - half_window):
            value = 0.0
            for offset, weight in offset_data:
                value += weight * paddata[i + offset][1]
            if value < 0:
                smooth_data.append((paddata[i][0], 0.0))
            else:
                smooth_data.append((paddata[i][0], value))
        return smooth_data

    def findpeaktops(self, data):
        """
        @brief finds all the local peak maxima in a spectrum
        @param data <list>: contains the MS spectrum data (mz <float>, inten <float>)
        @return maxima <list>: list of dictionaries containing the ion data of the maximum intensity data
        """

        pickParam = self.cfg.parameters['picking']
        if 'ms2valley' in pickParam and self.isMS2:
            valley = 1 - pickParam['ms2valley']
        else:
            valley = 1 - pickParam['valley']
        edge = pickParam['edge']
        thresh = 0
        maxima = []
        maxint = 0.0
        minint = 0.0
        last = 0
        for ion in xrange(len(data)):
            inten = data[ion][1]
            if inten > maxint:
                if maxima:
                    if maxima[-1]['right'] == 0:
                        maxima[-1]['right'] = ion - 1
                maxint = inten
                peak = dict(top=ion, topmz=data[ion][0], topint=inten, left=0, right=0, minleft=minint)
                thresh = peak['topint'] * edge
            elif inten < maxint * valley:
                if maxint > minint:
                    p = last

                    while data[p][1] < thresh:
                        p += 1
                    peak['left'] = p
                    maxima.append(peak)

                if inten < thresh and maxima[-1]['right'] == 0:
                    maxima[-1]['right'] = ion - 1
                maxint = inten
                minint = inten
                if maxima:
                    maxima[-1]['minright'] = inten
            elif inten <= minint:
                if maxint * valley > minint:
                    p = last
                    while data[p][1] < thresh:
                        p += 1
                    peak['left'] = p - 1
                    maxima.append(peak)

                if inten < thresh and maxima[-1]['right'] == 0:
                    maxima[-1]['right'] = ion - 1
                last = ion
                maxint = inten
                minint = inten
                if maxima:
                    maxima[-1]['minright'] = inten

        return maxima

    @staticmethod
    def centroid3(peak, spectrum):
        """
        @brief calculates the centroid of an ion using the top 3 data points
        @param peak <dictionary>: containing all the data about the peak
        @return mcentroid <float>: teh cantroid of the peak
        """
        left = spectrum[peak['top'] - 1]
        centre = spectrum[peak['top']]
        right = spectrum[peak['top'] + 1]
        mcentroid = (left[1] * left[0] + centre[1] * centre[0] + right[1] * right[0]) / (left[1] + centre[1] + right[1])
        return mcentroid

    @staticmethod
    def centroid(peak, spectrum):
        """
        @brief calculates the centroid of an ion using all data points in the peak
        @param peak <dictionary>: containing all the data about the peak
        @return mcentroid <float>: teh cantroid of the peak
        """
        sumint = 0.0
        sumproduct = 0.0
        for j in range(peak['left'], peak['right'] + 1):
            sumint += spectrum[j][1]
            sumproduct += spectrum[j][0] * spectrum[j][1]

        mcentroid = sumproduct / sumint
        return mcentroid

    @staticmethod
    def gaussian(peak, spectrum):
        """
        @brief calculates the gaussian fit of the peak using the top 3 data points
        @param peak <dictionary>: containing all the data about the peak
        @return mgauss <float>: the gaussian m/z of the peak
        """
        left = spectrum[peak['top'] - 1]
        centre = spectrum[peak['top']]
        right = spectrum[peak['top'] + 1]

        r2l = right[1] - left[1]
        c2r = centre[1] - right[1]
        l2c = left[1] - centre[1]

        mgauss = 0.5 * (c2r * left[0] ** 2 + r2l * centre[0] ** 2 + l2c * right[0] ** 2)
        mgauss /= (c2r * left[0] + r2l * centre[0] + l2c * right[0])
        return mgauss

    @staticmethod
    def loggaussian(peak, spectrum):
        """
        @brief calculates the gaussian fit of the peak using the top 3 data points and using log10 intensity data
        @param peak <dictionary>: containing all the data about the peak
        @return loggauss <float>: the log gaussian m/z of the peak
        """
        left = spectrum[peak['top'] - 1]
        centre = spectrum[peak['top']]
        right = spectrum[peak['top'] + 1]

        r = math.log10(right[1] + 1)
        l = math.log10(left[1] + 1)
        c = math.log10(centre[1] + 1)

        r2l = r - l
        c2r = c - r
        l2c = l - c
        loggauss = 0.5 * (c2r * left[0] ** 2 + r2l * centre[0] ** 2 + l2c * right[0] ** 2)
        loggauss /= (c2r * left[0] + r2l * centre[0] + l2c * right[0])
        return loggauss

    def filter4inten(self, peaks):
        """
        @brief filters the peaks by removing the lowest intensity fraction
        @param peaks <list>: list of all the peak data
        @return filtered peaks
        """
        filtered = peaks[:]
        filtered.sort(key=lambda x: x['maxinten'])

        num = int(len(filtered) * self.cfg.parameters['filtering']['fraction'])
        filtinten = filtered[num]['maxinten']
        filtered = filtered[num:]
        filtered.sort(key=lambda x: x['mz'])
        for f in range(len(filtered)):
            filtered[f]['index'] = f

        return filtered, filtinten

    def filter4charge(self, peaks):
        """
        @brief filters the peaks by removing all ions which can not be linked to any other ions as a potential isotope
        @param peaks <list>: list of all the peak data
        @return filtered peaks
        """
        mda = self.cfg.parameters['general']['tolmda']
        ppm = self.cfg.parameters['general']['tolppm']
        maxcharge = self.cfg.parameters['deisotoping']['max_charge']
        neutron = self.cfg.parameters['general']['neutron']
        filtered = []
        pending = {}

        for i in range(len(peaks)):
            base = peaks[i]
            error = max((mda, base['mz'] * ppm))
            offset = 0

            while 1:
                offset += 1
                if i + offset >= len(peaks):
                    break
                nextPeak = peaks[i + offset]
                delta = nextPeak['mz'] - base['mz']
                # stop if more than a neutron mass diff
                if delta > 1.1:
                    break

                charge = int(1 / delta + 0.5)
                thmz = base['mz'] + neutron / charge
                if abs(nextPeak['mz'] - thmz) > error or charge > maxcharge:
                    # doesnt fit to charge so skip to the next ion
                    continue
                pending[i] = 1
                pending[i + offset] = 1

            if i in pending:
                idx = len(filtered)
                filtered.append(peaks[i].copy())
                filtered[-1]['index'] = idx
                del pending[i]

        if pending:
            for k in pending:
                filtered.append(peaks[k].copy())

        return filtered


class IsotopeModel:
    """
    @brief class to manage theoretical isotope distributions of peptides
    """

    def __init__(self, cfg, logs, maxrepeats, mininten, quanmeth, overwrite=0):
        """
        @brief initiates the averagine isotopic abundance model
        @param cfg <object>: contains all the configuration data
        @param maxrepeats <integer>: maximum number of averagines to calculate the isotopic abundances for
        @param mininten <float>: minimum fraction of the normalised intensity for an isotope to be used
        @param quanmeth <string>: the quantitation method used to select the correct isofile
        """
        self.deisoparam = cfg.parameters['deisotoping']
        self.maxRepeats = maxrepeats
        self.minInten = mininten
        self.c12 = 0
        self.isolist = None
        self.numlabels = self.deisoparam['num_labels']
        self.heavyIsotopes = ['C13', 'N15', 'O18', 'H2']
        if quanmeth:
            methName = quanmeth['meth_name']
            self.isotopeNumbers = quanmeth['heavy_isotopes']
        else:
            methName = 'none'
            self.isotopeNumbers = {}

        isof = Path('./data/averagine %s.txt' % methName.lower())

        if not isof.parent.exists():
            isof.parent.mkdir()
        if not isof.exists() or overwrite:
            self.calcIsotopeDistribution(methName)
            self.savedata(isof, methName)
            logs.log.info('created isotope data and saved to "%s"' % isof.name)
        else:
            self.loadfile(isof)
            logs.log.info('loaded isotope data from "%s"' % isof.name)
            # self.calcIsotopeDistribution()
            # self.savedata()

    def loadfile(self, isof):
        """
        @brief loads the isotope ratios from file
        @param isof <object>: path object for the isoratio file
        """
        f = file(str(isof), 'r')
        line = f.next()
        splt = line[:-1].split()
        limits = (float(splt[0]), int(splt[1]))
        self.c12 = int(splt[2])
        self.maxRepeats = limits[1]
        isos = []
        for j in range(limits[1]):
            line = f.next()
            data = line[:-1].split('\t')
            mass = float(data[0])
            comp = data[1]
            ratios = [float(x) for x in data[2:]]
            iso = dict(monomass=mass, comp=comp, inten=ratios)
            isos.append(iso)
        self.isolist = isos

    def findIsotopeSet(self, mass, remedge=0):
        """
        @brief finds the closest averagine model to the mass
        @param mass <float>: the mass to find the averagine for
        @param remedge <integer>: flag to remove isotopes below the 12C, 1 removes below 12C 2 also removes high
        isos according to mass
        @return iso <tuple>: the averagine isotope model
        (m/z <float>, num isos <integer>, isotope ratios <list of floats>)
        """
        a = 0
        while a < self.maxRepeats:
            if mass < self.isolist[a]['monomass']:
                break
            a += 1
        if a == self.maxRepeats:
            iso = self.isolist[-1].copy()
        elif abs(mass - self.isolist[a]['monomass']) <= abs(mass - self.isolist[a - 1]['monomass']):
            iso = self.isolist[a].copy()
        else:
            iso = self.isolist[a - 1].copy()

        # find the isotope with the greatest intensity
        maxpos = 0
        maxint = 0
        for i in range(len(iso['inten'])):
            if iso['inten'][i] > maxint:
                maxint = iso['inten'][i]
                maxpos = i
        iso['max'] = maxpos - self.c12

        if remedge >= 1:
            iso['inten'] = iso['inten'][self.c12:]
            maxpos -= self.c12

        if remedge == 2:
            edge = maxpos + self.deisoparam['num_past_max'] + 1
            if edge > 4:
                iso['inten'] = iso['inten'][:edge]
            else:
                iso['inten'] = iso['inten'][:4]

        return iso

    def calcIsotopeDistribution(self, quanmeth, maxrepeats=0, mininten=0):
        """
        @brief calculates the averagine isotopic abundances
        @param quanmeth <string>: the quantitation method used to select the correct isofile
        """
        # using the average of 20 amino acids
        # elem = dict(C=(5.35, 12.0, 12.0107,  0.011122346),
        #             H=(7.85, 1.007825032, 1.00794, 0.00015),
        #             N=(1.45, 14.0030740052, 14.0076, 0.003713741),
        #             O=(1.45, 15.9949146221, 15.9994, 0.002, 0.000380914),
        #             S=(0.10, 31.97207069, 32.065, 0.044306462, 0.007893075))
        # using the average of 20 amino acids weighted by database frequency
        elem = dict(C=(4.91, 12.0, 12.0107, 0.011122346),
                    H=(7.75, 1.007825032, 1.00794, 0.00015),
                    N=(1.38, 14.0030740052, 14.0076, 0.003713741),
                    O=(1.48, 15.9949146221, 15.9994, 0.000380914, 0.002),
                    S=(0.04, 31.97207069, 32.065, 0.007893075, 0.044306462)
                    )
        if maxrepeats:
            self.maxRepeats = maxrepeats
        if mininten:
            self.minInten = mininten

        self.isolist = []
        iso = dict(monomass=0.0, avgmass=0.0, inten=[], norm=[], C=0, H=0, N=0, O=0, S=0, C13=0, N15=0, O18=0, H2=0)

        # start with the labels if any
        if quanmeth == 'none':
            sumiso = [1]
            self.c12 = 0
        else:
            sumiso, self.c12 = self.labeldist(self.isotopeNumbers, self.numlabels, iso)

        for rep in range(1, self.maxRepeats + 1):
            # add the isotope data
            for e in elem:
                for atom in range(int(iso[e] + 0.5), int(iso[e] + elem[e][0] + 0.5)):
                    iso['monomass'] += elem[e][1]
                    iso['avgmass'] += elem[e][2]

                    sumiso.append(0)
                    if e in 'OS':
                        sumiso.append(0)

                    for i in range(len(sumiso) - 1, 0, -1):
                        sumiso[i] += sumiso[i - 1] * elem[e][3]

                        if e in 'OS' and i > 1:
                            sumiso[i] += sumiso[i - 2] * elem[e][4]

                iso[e] += elem[e][0]

            # normalise the data to the maximum intensity
            maxint = max(sumiso)
            c12int = sumiso[self.c12]
            iso['norm'] = []
            iso['inten'] = []
            maxnomalsedint = 0
            for i in sumiso:
                normint = i / maxint
                if normint > maxnomalsedint:
                    maxnomalsedint = normint

                if normint >= self.minInten:
                    iso['norm'].append(normint)
                    iso['inten'].append(i / c12int)
                elif maxnomalsedint < 1:
                    iso['norm'].append(0.0)
                    iso['inten'].append(i / c12int)

            self.isolist.append(iso.copy())
        return 1

    def savedata(self, isof, quanmeth):
        """
        @brief saves the calculated isotope data to file
        @param isof <object>: path object for the isoratio file
        @param quanmeth <string>: the quantitation method used to select the correct isofile
        """
        fin = open(str(isof), 'w')

        fin.write('%.3f\t%d\t%d\t%s\t%d\n' % (self.minInten, self.maxRepeats, self.c12, quanmeth, self.numlabels))

        for iso in self.isolist:
            line = '%.6f' % iso['monomass']
            line += '\tC%d H%d N%d O%d S%d' % (int(iso['C'] + 0.5), int(iso['H'] + 0.5),
                                               int(iso['N'] + 0.5), int(iso['O'] + 0.5), int(iso['S'] + 0.5))

            for label in self.heavyIsotopes:
                if iso[label] > 0:
                    line += ' %s %d' % (label, iso[label])
            for i in iso['inten']:
                line += '\t%.4f' % i

            fin.write(line + '\n')
        fin.close()

    def labeldist(self, label, rpts, iso):
        """
        @brief calculates the label isotope distribution
        @param label <string>: the type of isotopic label used in the experiment
        @param rpts <integer>: the maximum number of averagines used in the calculations
        @param iso <dictionary>:  the isotope dictionary
        @return niso , c12
        """
        isoData = dict(C13=[13.0033548378, 12.9933, 0.008],
                       N15=[15.0001088984, 14.9901, 0.008],
                       O18=[17.9991610700, 17.9833, 0.008],
                       H2=[2.014101778, 2.0061, 0.008])
        elem = {}
        for e in isoData:
            if e in label:
                elem[e] = [label[e] * rpts] + isoData[e]
            else:
                elem[e] = [0] + isoData[e]

        sumiso = [1]

        # add the labeled isotope data
        for e in self.heavyIsotopes:
            num = int(elem[e][0] + 0.5)
            for atom in range(1, num + 1):
                iso['monomass'] += elem[e][1]
                iso['avgmass'] += elem[e][2]

                sumiso = [0] + sumiso

                for i in range(len(sumiso) - 1):
                    sumiso[i] += sumiso[i + 1] * elem[e][3]

            iso[e] += num

        # filter the isotopes
        div = max(sumiso)
        niso = []
        for i in sumiso:
            norm = i / div
            if norm == 1:
                c12 = len(niso)

            if norm > 0.001:
                niso.append(norm)

        return niso, c12
