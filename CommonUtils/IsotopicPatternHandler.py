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

This file mainly deals with the handling of the command line data and reading
configuration files converting parameter values as directed.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob

"""

import numpy as np

class IsotopicPatternHandler:
    def __init__(self, hdf5Object, cfg, XICprocessor):
        self.hdf5 = hdf5Object
        self.cfg = cfg
        self.XICprocessor = XICprocessor

    def getMeasuredClustersFromXIC(self, theoretical, pcm):
        cfg = self.cfg

        error = max((cfg.parameters['deisotoping']['tolmda'],
                     theoretical['MZs'][0] * cfg.parameters['deisotoping']['tolppm']))
        minRT = pcm.RTmsms - cfg.parameters['xic']['afterpeak']
        maxRT = pcm.RTmsms + cfg.parameters['xic']['beforepeak']
        allPeaks = {}

        usedRelInts = {}
        if cfg.parameters['deisotoping']['limit_xics']:
            # if limitXICs:
            for i in range(-cfg.parameters['deisotoping']['nummissincorp'],
                           cfg.parameters['deisotoping']['numpeaksforfit']):
                usedRelInts[i] = theoretical['primary'][i]
        else:
            usedRelInts = theoretical['primary']

        isoIndexes = usedRelInts.keys()
        lowMZindex = min(isoIndexes, -1)
        highMZindex = max(isoIndexes)

        # load small XIC block for all isotopes
        self.XICprocessor.setXICblock(theoretical['MZs'][lowMZindex], theoretical['MZs'][highMZindex], minRT, maxRT)

        # generate XIC peaks for all isotope channels
        for idx in usedRelInts:
            mz = theoretical['MZs'][idx]

            xicPts = self.XICprocessor.createXIC(mz, error)

            allPeaks[idx] = []  # commented for debug reasons
            if xicPts:
                peaks = self.XICprocessor.findPeaksInXIC(minRT, maxRT, 0)
                ##speaks = self.XICprocessor.findPeaksInXICNew(minRT, maxRT, 0)

                # print '%i peaks for [%i]: %.6f  (%.1f-%.1f)' % (len(peaks), idx, mz, minRT, maxRT)
                for p in peaks:
                    # if p[0]['rt'] >= minRT and p[0]['rt'] <= maxRT:
                    allPeaks[idx].append(p[0])

        # find prior ion data
        priorIonMZ = theoretical['MZs'][-1]
        xicPts = self.XICprocessor.createXIC(priorIonMZ, error)

        # fout = open('c:/fs/data/priorIonData.txt', 'a')

        for peak in allPeaks[0]:
            # use all the 12C peaks and find the prior ion data
            priorinten, apexspec = self.XICprocessor.findMaxIntInTimeRange(peak['rt50left'], peak['rt50right'], peak['rt'])
            #self.
            #priorInten = self.XICprocessor.findMaxIntInTimeRangeOLD(peak['rt50left'], peak['rt50right'], peak['rt'])
            peak['priorinten'] = priorinten
            peak['apexSpec'] = apexspec
            #self.XICprocessor.apexSpec

            # export data for output
          #  peak['datalist'] = self.XICprocessor.dataList
        #     if self.XICprocessor.dataList:
        #
        #         for point in self.XICprocessor.dataList:
        #             fout.write('%i\t%f\t%f\t%f\t%i\t%f\t%f\n' % (peak['specid'], peak['peakmz'], peak['rt'],
        #                        peak['inten'], point[0], point[1], point[2]))
        #     else:
        #         fout.write('%i\t%f\t%f\t%f\t0\t0\t0\n' % (peak['specid'], peak['peakmz'], peak['rt'], peak['inten']))
        #
        #     x = 1
        # fout.close()

        clusters = []
        for monoPK in allPeaks[0]:
            # iterate through the peaks in the mono-isotopic peak channel to find more cluster data
            cluster = {0: monoPK.copy()}
            for iso in allPeaks:
                if iso == 0:
                    continue
                for isoPK in allPeaks[iso]:
                    # compare RT with monoisotopic data
                    if monoPK['rt50left'] <= isoPK['rt50right'] and monoPK['rt50right'] >= isoPK['rt50left']:
                        # peak matches
                        if iso in cluster:
                            # already a fitted peak for that isotope
                            if abs(isoPK['rt'] - monoPK['rt']) < abs(cluster[iso]['rt'] - monoPK['rt']):
                                cluster[iso] = isoPK.copy()
                        else:
                            cluster[iso] = isoPK.copy()

            # all peaks tried
            cluster['origin'] = 'XIC'
            # filter data for the presence of 12C and 13C data
            if 0 in cluster and 1 in cluster:
                cluster['priorinten'] = cluster[0]['priorinten']
                clusters.append(cluster.copy())

        return clusters

    def getMeasuredClustersFromSpectra(self, spectrum, theoretical):
        """
        @param spectrum:
        @param theoretical:
        @return: found - a ghastly dictionary with some values and some nested vals.  should just be an object.
        """
        clusters = []
        error = max(theoretical['monoMZ'] * self.cfg.parameters['deisotoping']['tolppm'],
                    self.cfg.parameters['deisotoping']['tolmda'])
        largeError = max(theoretical['monoMZ'] * 50e-6, 0.05)
        # code below yields a matrix of len(spectrum) with all FALSE values (we expect
        # no MZs with zero values! This is the used to record where we find MZ matches within tolerance to the
        # theoretical values as given in theoretical
        foundmatrix = spectrum['mz'] == 000
        found = dict(MZs={})
        found_isos = 0
        missing = dict()
        for i, theomz in theoretical['MZs'].items():
            matchmatrix = abs(spectrum['mz'] - theomz) < error
            if matchmatrix.any():  # match within tolerance
                best = sorted(spectrum[matchmatrix], key=lambda x:abs(x[3]-theomz))[0]
                found['MZs'][i] = best['mz']
                maxinten = best['inten']
                found[i] = maxinten
                if i == 0:
                    found['specMZ'] = best['mz']
                    found['specRT'] = best['rt']
                    found_isos += 1
                elif i == 1:
                    found_isos += 1
                foundmatrix = matchmatrix | foundmatrix # add the position of the matched mzs to the large matrix
            else:
                # to record the missing data (c12 / c13 only for interest's sake) get the abs delta to the next value
                # if within largeError tolerance
                if i in (0, 1):
                    missingmatchmatrix = abs(spectrum['mz'] - theomz) < largeError
                    if missingmatchmatrix.any():
                        missing[i] = (spectrum[missingmatchmatrix]['mz'][0], spectrum[missingmatchmatrix]['inten'][0],
                                      spectrum[missingmatchmatrix]['mz'][0] - theomz)

        found['origin'] = "spec"
        if found_isos == 2:
            clusters.append(found)
            return clusters, {}
        else:
            return clusters, missing

    @staticmethod
    def correctForOverlappingIsoPatterns(clusters, allTheoretical, label, prevLabel, clusterMatches, overlapInfluence):
        """

        @param clusters: warning cluster is changed in place and also returned!
        @param allTheoretical:
        @param label:
        @param prevLabel:
        @param clusterMatches:
        @return:
        """
        # skip = ['origin', 'priorinten', 'usedInten', 'ovelapinfluence']
        #print 'starting correct for overlapping', label, prevLabel
        theoretical = allTheoretical[label]

        for cluster in clusters:
            # extract intensity data
            clustInten = {}

            if cluster['origin'] is "XIC":
                isos = [x for x in cluster if isinstance(cluster[x], dict)]
                isos.sort()

                clustInten[-1] = cluster['priorinten']
                for iso in isos:
                    clustInten[iso] = cluster[iso]['smooth']
                isos.insert(0, -1)
            elif cluster['origin'] is "spec":
                isos = [x for x in cluster if isinstance(x, int)]
                isos.sort()
                for iso in isos:
                    clustInten[iso] = cluster[iso]
            else:
                raise ValueError((cluster['origin'] + ' is not a valid value for cluster[\'origin\']'))

            if prevLabel and prevLabel in clusterMatches:
                overlapsFromLighter = theoretical['overlaps'].count('2')

                if overlapsFromLighter > 0:
                    prevTheo = allTheoretical[prevLabel]
                    prevMatched = clusterMatches[prevLabel]

                    endInCurrent = min(theoretical['primary'].keys()) + overlapsFromLighter
                    peakInPrev = max(prevTheo['primary'].keys()) - (overlapsFromLighter - 1)

                    for peakInCurrent in range(min(theoretical['primary'].keys()), endInCurrent):

                        if peakInCurrent in cluster or (peakInCurrent == -1 and cluster['origin'] is "XIC"):
                            tmpInten = clustInten[peakInCurrent]
                            if 'fitInten' in prevMatched:
                                tmpInten -= prevMatched['fitInten'] * prevTheo['primary'][peakInPrev] * overlapInfluence

                            if tmpInten > 0:
                                clustInten[peakInCurrent] = tmpInten
                            else:
                                clustInten[peakInCurrent] = 1
                            pass

                        peakInPrev += 1

            cluster['usedInten'] = clustInten.copy()

        return clusters

    def getBestFittingIsoPatternCandidate(self, isoPatternCandidates, theoretical, charge, useSecondary):
        """

        @param isoPatternCandidates:
        @param theoretical:
        @param averagine:
        @return:
        """

        bestFit = 10000
        bestCluster = {}
        p1 = 0
        for patternCand in isoPatternCandidates:
            primaryWanted = self.getPeaksForMatching(theoretical['primary'])
            if useSecondary:
                secondaryWanted = self.getPeaksForMatching(theoretical['secondary'])

            if 'usedInten' not in patternCand:
                continue

            clustInten = patternCand['usedInten']

            p1 += 1
            maxInten = max(clustInten.values())
            if maxInten > 1 and clustInten[0] > 1:
                self.cfg.modeltype = 'ext'
                fit = self.isofitleastsq(primaryWanted, clustInten, charge, 1)
                if useSecondary:
                    self.cfg.modeltype = 'avg'
                    secondaryFit = self.isofitleastsq(secondaryWanted, clustInten, charge, 1)
                else:
                    secondaryFit = [None, None, None]

                if fit[0] < bestFit:
                    bestFit = fit[0]
                    bestCluster = dict(cluster=patternCand.copy(), leastSquares=fit[0], isoMatched=fit[1],
                                       fitInten=fit[2], secondary_fit=secondaryFit[0], secondary_inten=secondaryFit[2])

        if not isoPatternCandidates:
            pass

        return bestCluster

    def getPeaksForMatching(self, theoretical):
        """
        @brief filters the theoretical isotopes to only those needed for matching
        @param theoretical <dictionary>:  containing the thoeretical intensities indexed by the 13C isotope numbers
        @return filtered <dictionary>: the filtered intensity dictionary
        """

        filtered = {}
        for i in range(-self.cfg.parameters['deisotoping']['nummissincorp'],
                       self.cfg.parameters['deisotoping']['numpeaksforfit']):
            if i in theoretical:
                filtered[i] = theoretical[i]

        return filtered

    def isofitleastsq(self, theoretical, measured, charge, havestart):
        """
        @brief find the best isotope fit to the ion data using least squares optimisation
        @param theoretical <list>: of the relative intensities of the theoretical isotope model
        @param measured <dictionary>: containing the observed ion data indexed by the isotope offset
        @param charge <integer>: charge state of the ion to be deisotoped
        @param havestart <integer>: flag indicating if the data should have the start ion present
        @return ion <tuple>: containing the mz, intensity and offset of the 12C ion of the preferred fit
        """
        calcsumsq = self.calcsumsq
        sellowest = self.sellowest
        # cfg = self.cfg
        neutron = self.cfg.parameters['general']['neutron']

        avint = 0.0
        num = len(theoretical)

        fits = []
        intensities = {}
        minint = 1e10
        maxint = 0
        matchdic = {}
        summass = 0.0
        summassint = 0.0
        sumint = 0.0
        # ensure that the list of theoretical masses is the same as measured.
        for i in theoretical:
            if i in measured:
                inten = measured[i]
                matchdic[i] = 1
            # summass += (measured[start + i][0] - i * neutron / charge)
            # summassint += (measured[start + i][0] - i * neutron / charge) * measured[start + i][1]
            # sumint += measured[start + i][1]
            else:
                inten = 1.0
                matchdic[i] = 0
            if inten > maxint:
                maxint = inten
            if inten < minint:
                minint = inten
            intensities[i] = inten

        if maxint > 0:
            opts = []
            next = -1
            opts.append((minint, calcsumsq(theoretical, intensities, minint)))
            midint = (maxint + minint) / 2
            opts.append((midint, calcsumsq(theoretical, intensities, midint)))
            opts.append((maxint, calcsumsq(theoretical, intensities, maxint)))
            if opts[0][1] < opts[1][1]:
                # first point is lowest - need another point lower
                newint = opts[0][0] / 2
                opts.append((newint, calcsumsq(theoretical, intensities, newint)))
                opts.sort()
                opts, next = sellowest(opts)

            if opts[2][1] < opts[1][1]:
                # last point is lowest - need another point
                newint = opts[2][0] * 2
                opts.append((newint, calcsumsq(theoretical, intensities, newint)))
                opts.sort()
                opts, next = sellowest(opts)

            if opts[0][1] < opts[1][1] or opts[2][1] < opts[1][1]:
                # single extension of range has not worked - cancel fitting
                return 10000, {0: 0}, 0, 0, 0

            if next == -1:
                opts, next = sellowest(opts)

            start = opts[:]
            loop = 0
            while (opts[next][1] - opts[1][1]) / opts[1][1] > self.cfg.parameters['deisotoping']['ls_limit'] and \
                  (opts[2][0] - opts[0][0]) / opts[1][0] > self.cfg.parameters['deisotoping']['int_limit']:
                midint = (opts[1][0] + opts[next][0]) / 2
                opts.append((midint, self.calcsumsq(theoretical, intensities, midint)))
                opts.sort()
                opts, next = sellowest(opts)
                loop += 1
                # fits.append((opts[1][1] / opts[1][0], matches, opts[1][0], start))
            fits.append((opts[1][1], matchdic, opts[1][0]))
        if fits[0][0]<0.1:
            #print 'measured', [x / float(measured[0]) for x in measured.values()]
            #print 'theoretical', [x for x in theoretical.values()]

            lastval = min(len(measured.values()),len(theoretical.values()))
           # print 'm-t',self.cfg.mypcm, self.cfg.thislab,self.cfg.model,self.cfg.modeltype, np.array([x / float(measured[0]) for x in measured.values()[:lastval]]) - np.array( [x for x in theoretical.values()][:lastval]),
           # print 'fits',fits[0][0]

        if fits:
            fits.sort()
            return fits[0]
        else:
            return 10000, {0: 0}, 0, 0, 0

    @staticmethod
    def calcsumsq(theoretical, measured, fitIntensity):
        sumsq = 0.0
        for i in theoretical:
            sumsq += (theoretical[i] - measured[i] / fitIntensity) ** 2
        return sumsq

    @staticmethod
    def sellowest(data):
        t = data[:]
        if len(data) < 4:
            pass
        elif data[1][1] < data[2][1]:
            data = data[:3]
        else:
            data = data[1:]

        if data[0][1] > data[2][1]:
            next = 0
        else:
            next = 2
        return data, next

    @staticmethod
    def getMisincorporationRatios(associatedIsoPatterns, theoretical):
        for label in theoretical['byIncreasingMZ']:

            xicForLabel = associatedIsoPatterns['XIC'][label]
            apexForLabel = associatedIsoPatterns['apexSurvey'][label]
            priorForLabel = associatedIsoPatterns['priorSurvey'][label]
            misincorporationRatio_apex = 0
            misincorporationRatio_prior = 0

            if 'cluster' in xicForLabel:
                if 'priorinten' in xicForLabel['cluster']:
                    xicCluster = xicForLabel['cluster']
                    misincorporationRatio_old = xicCluster['priorinten'] / (theoretical[label]['primary_sumTheo'] *
                                                                            xicForLabel['fitInten'])
                    if xicCluster['usedInten'][-1] == 1:
                        xicCluster['usedInten'][-1] = 0
                    misincorporationRatio_xic = xicCluster['usedInten'][-1] / (theoretical[label]['primary_sumTheo'] *
                                                                               xicForLabel['fitInten'])
            if 'cluster' in apexForLabel:
                if -1 in apexForLabel['cluster']:
                    apexCluster = apexForLabel['cluster']
                    misincorporationRatio_apex = apexCluster['usedInten'][-1] / (theoretical[label]['primary_sumTheo'] *
                                                                                 apexForLabel['fitInten'])

            if 'cluster' in priorForLabel:
                if -1 in priorForLabel['cluster']:
                    priorCluster = priorForLabel['cluster']
                    misincorporationRatio_prior = priorCluster['usedInten'][-1] / (theoretical[label]['primary_sumTheo']
                                                                                   * priorForLabel['fitInten'])

            if associatedIsoPatterns['XIC'][label]:
                associatedIsoPatterns['XIC'][label]['-1ratio'] = misincorporationRatio_xic
                associatedIsoPatterns['XIC'][label]['-1ratio_old'] = misincorporationRatio_old

            if associatedIsoPatterns['apexSurvey'][label]:
                associatedIsoPatterns['apexSurvey'][label]['-1ratio'] = misincorporationRatio_apex

            if associatedIsoPatterns['priorSurvey'][label]:
                associatedIsoPatterns['priorSurvey'][label]['-1ratio'] = misincorporationRatio_prior
