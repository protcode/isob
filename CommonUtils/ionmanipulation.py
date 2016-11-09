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

This file mainly performs MS data manipulation e.g. deconvolution or deisotoping

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""


class _Cfg:
    """config container"""
    def __init__(self):
        pass


class IonProcessing:
    """class for methods used for processing MS2 ion data
    """

    def __init__(self, cfg, log):

        self.cfg = cfg
        self.log = log
        self.masses = dict(Proton=1.0072760000000001)

    def doIsotopeCorrection(self, correctionFactors, valDict):
        """
        @brief For each reporter ion, a proportion of the total signal of the other reporter ions make
        it into the signal due to differing amounts of C13.  These proportions should be subtracted from each
        signal
        @param correctionFactors <dict> For all other colours given for the key, the proportion
        of the colour's total signal which making up the total of key's signal is given
        @param validict <dict>: The reporter ion signals to be corrected given by each key
        @return correctedData
        """
        # correctvals is a dictionary of matries per colour
        self.log.debug('started doIsotopeCorrection with %s, %s' % (str(correctionFactors), str(valDict.keys())))
        interferenceData = {}

        for iso in valDict:
            # calculate the total interference for each potential ion
            self.log.debug('extracting interference values for colour:  %i' % iso)
            interference = 0
            for id, corrFactor in correctionFactors[iso].items():
                interference = valDict[iso]['area'] * corrFactor
                self.log.debug('proportion of signal from %i to %i (ratio %s of total %s) is %s '
                               % (iso, id, corrFactor, valDict[iso]['area'], interference))
                if interference < 0:
                    self.log.debug('interference below zero: setting to zero')
                    interference = 0

                try:
                    interferenceData[id] += interference
                except KeyError:
                    interferenceData[id] = interference

        for iso in valDict:
            # apply the interference values to the detected isotopes
            try:
                corrected = valDict[iso]['area'] - interferenceData[iso]
            except KeyError:
                corrected = valDict[iso]['area']

            if corrected < 0:
                valDict[iso]['corrected'] = 0
            else:
                valDict[iso]['corrected'] = corrected

        self.log.debug('interferenceData %s ' % str(interferenceData))
        return

    def deisotope(self, header, ions):
        """
        @brief removes 13C isotopes from the MS/MS spectrum
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        @return removed <integer>: the number of ions removed from the spectrum
        """
        neutron = self.cfg.parameters['general']['neutron']

        removed = 0
        # calculate the possible deltas for the precursor charge
        deisotopeda = 0.02

        if header['fragmeth'].upper() == 'HCD':
            neutrondelta = [(neutron / c, c) for c in range(header['charge'], 0, -1)]

        else:
            # all other ms methods should not be so precise and cannot guarantee
            # different charge states, so just check for 1+
            neutrondelta = [(neutron, 1)]

        # loop backwards through the data
        for c13 in range(len(ions) - 1, 0, -1):
            self.log.debug('c13 %s, neutrondelta %s' % (c13, neutrondelta))
            mz13c = ions[c13]['mz']
            if mz13c < 132:
                break
            c12 = c13
            delta = 0
            # now scan back for a prior ion
            while delta < 2 and c12 > 0:  # no point looking if the difference over 2or if the C12 index reached zero

                c12 -= 1
                delta = mz13c - ions[c12]['mz']
                self.log.debug('for spectrum %s : delta %s,  ions[c12] %s ions[c13] %s' %
                               (header['spec_id'], delta, ions[c12]['mz'], ions[c13]['mz']))
                for nd, cge in neutrondelta:
                    # if delta minus neutron within tolerance (deisotopeda and c12 intensity more than half of c13
                    # can remove c13 peak and calculate charge of fragment
                    if abs(delta - nd) < deisotopeda and ions[c12]['inten'] > (ions[c13]['inten'] / 2.0):
                        self.log.debug('allowing deisotoping %s deisotopeda %s delta %s c12 %s c13 %s' %
                                       (nd, deisotopeda, delta, ions[c12]['mz'], ions[c13]['mz']))
                        # ion matches in mass and intensity
                        ions[c12]['charge'] = str(cge)
                        ions[c13]['charge'] = str(cge)
                        ions[c13]['inten'] = 0
                        removed += 1

        return ions, removed

    def deconvolute(self, header, ions):
        """
        @brief deconvolutes the MS/MS spectrum
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        """
        proton = self.masses['Proton']

        # remove 13C ions and identify charge states
        ions, removed = self.deisotope(header, ions)
        if removed:
            # only do the deconvolution if data was removed ie isotopes were detected
            for ion in ions:
                cge = float(ion['charge'])
                self.log.debug('cge %s ion %s' % (cge, ion))
                # if we found charge for ion we can deconvolute it
                if cge > 1:
                    deconv = (ion['mz'] - proton) * cge + proton
                    self.log.debug('header %s mz %s deconv %s' % (header['spec_id'], ion['mz'], deconv))
                    parentmass = header['precmz'] * header['charge']
                    self.log.debug('parentmass %s' % parentmass)
                    if deconv > parentmass > 0:
                        ion['inten'] = 0
                    else:
                        ion['mz'] = deconv
        ions.sort(order='mz')
        # 1.5 is abitrary value to ensure the missincorp of label is never seen as an ion
        maxmass = (header['monomz'] - proton) * header['charge'] - 1.50
        idx = ions['mz'] < maxmass
        return ions[idx]
