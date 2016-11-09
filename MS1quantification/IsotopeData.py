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
# python imports
from operator import itemgetter

# CommonUtils imports

# project imports


class IsotopeData:
    """
    @brief storage of elemental compositions and isotope patterns for an elemental composition and its modifications for
            ms1 quantification
    """

    protonMass = 1.00727645
    neutronMass = 1.0033548378
    minPeakDist = 35

    fixedModsIsoDist = None
    ms1QuantStates = []


    def __init__(self, elemComp, sequence, ms1QuantMods, ECC, IC, MS1QuantOp, exactUseage,
                 isotopeCalculationMethod='bySEQ', varMods=None):
        """
        @brief initialize with a light elemental composition, the corresponding peptide sequence
                and one instance of each ElementalCompositionCalculator() and IsotopeCalculator()
                elemental compositions and isotope patterns modified for ms1 quantification are automatically retrieved
                and stored
        @param elemComp: <dictionary> light elemental composition
        @param sequence: <string> corresponding peptide sequence
        @param ECC: instance of ElementalCompositionCalculator()
        @param IC: instance of IsotopeCalculator()
        @param MS1QuantOp: instance of MS1QuantOperations()
        """

        # set the instance variables
        self.sequence = sequence
        self.icMethod = isotopeCalculationMethod
        self.massesAndRelInts = None
        self.elemComp = elemComp
        self.massesAndRelIntsByCharge = {}
        self.exactUsage = exactUseage

        # if the ms1QuantStates class variable is not yet set, add the appropriate states to it
        if not IsotopeData.ms1QuantStates:
            for sm in ms1QuantMods:
                if not sm['label'] in IsotopeData.ms1QuantStates:
                    IsotopeData.ms1QuantStates.append(sm['label'])

        # if isotopic distributions are to be computed based on elemental compositions, generate modified elemental
        #   compositions
        if self.icMethod is 'byEC_noClip' or self.icMethod is 'byEC_clip':
            modifiedECs = MS1QuantOp.applyMS1QuantModToElementalComposition(sequence, elemComp)

        # call the slow version of elemental composition-based calculation of isotopic distributions
        if self.icMethod == 'byEC_noClip':
            self.MS1QuantIsoPatterns = {}

            # loop through the labelling sates
            for state in IsotopeData.ms1QuantStates:
                # for the modified versions, call the composition-based calculation of isotopic distributions with the
                #   appropriately modified elemental compositions
                self.MS1QuantIsoPatterns[state] = IC.filterRelIntensities(IC.calcIsotopeDistributionFromElemComp(
                    modifiedECs[state]['elemComposition']), False)

        # call the fast version of elemental composition-based calculation of isotopic distributions (still slower than
        #   sequence-based)
        elif self.icMethod == 'byEC_clip':
            self.MS1QuantIsoPatterns = {}

            # loop through the ms1 quant sates
            for state in IsotopeData.ms1QuantStates:
                # for the modified versions, call the composition-based calculation of isotopic distributions with the
                #   appropriately modified elemental compositions
                self.MS1QuantIsoPatterns[state] = IC.filterRelIntensities(IC.calcIsotopeDistributionFromElemComp(
                    modifiedECs[state]['elemComposition']))

        # call the sequence-based version for calculation of isotopic distributions (fastest version)
        elif self.icMethod == 'bySEQ':
            # strip out the modified amino acids from the original sequence and calculate the isotopic distribution
            #   of the resulting 'core sequence'
            coreSequence = MS1QuantOp.getCoreSequence(sequence)
            coreIsoPattern = IC.calcIsotopeDistributionFromSequence(coreSequence)

            # if the are variable modification defined for the current peptide, apply them
            if sum(varMods.values()):
                # calculate the isotopic distribution of the pure variable modifications and combine it with that of
                #   the core sequence
                varModsIsoPattern = IC.calcIsotopeDistributionFromElemComp(varMods)
                coreIsoPattern = IC.addIsotopePatterns(coreIsoPattern, varModsIsoPattern)

            # prepare a dictionary for ms1 labeled isotope patterns
            self.MS1QuantIsoPatterns = {}

            # loop through the ms1 quant sates
            for state in IsotopeData.ms1QuantStates:
                # for the modified versions, get the isotope pattern of the correspondingly modified amino acids
                #   and combine it with that of the core sequence
                complementPattern = MS1QuantOp.getComplementaryIsopattern(sequence, state)
                self.MS1QuantIsoPatterns[state] = IC.filterRelIntensities(
                    IC.addIsotopePatterns(coreIsoPattern, complementPattern))

    def ppmDistance(self, mzA, mzB):
        """
        @brief compute the difference between two mz values in ppm
        @param mzA: <float> one of the mz values
        @param mzB: <float> the other mz value
        @return: difference between mzA and mzB in ppm (relative to mzA)
        """
        return (mzA - mzB) / mzA * 1000000

    def getMassesAndRelInts(self, charge):
        """
        @brief For a given charge, return the mz values and corresponding relative intensities for the isotopic patterns
                of all ms1-quant modified states of a PM
        @param charge: <int> the charge needed to compute the mz values
        @return: dictionary containing the mz values and relative intensities of the isotopic peaks for each labelling
                   state according to the given mass
        """
        neutron = IsotopeData.neutronMass
        relInts = self.exactUsage
        sumTheo = relInts + '_sumTheo'

        # if the results for the given charge have not yet been computed
        if charge not in self.massesAndRelIntsByCharge:
            # prepare a new dictionary for the results
            massesAndRelInts = {}

            # loop through the labelling states
            for state in self.MS1QuantIsoPatterns:
                c12Idx = self.MS1QuantIsoPatterns[state]['c12index']
                # create a entry in the return dictionary for monoisotopic mz and compute it
                monoMZ = (self.MS1QuantIsoPatterns[state]['monomass'] + (charge * IsotopeData.protonMass)) / charge
                massesAndRelInts[state] = {'monoMZ': monoMZ}

                # prepare lists for mz values, relative intensities and detected overlaps between the ms1-quant isotope
                #   patterns
                massesAndRelInts[state]['MZs'] = {}
                massesAndRelInts[state][relInts] = {}
                massesAndRelInts[state]['overlaps'] = []

                # loop through the peaks in the isotope pattern of the current ms1-quant state
                for peakNo in xrange(len(self.MS1QuantIsoPatterns[state]['norm'])):
                    # prepare the overlaps list by appending a '0' for each peak
                    # (indicating no overlap, may be changed later on)
                    massesAndRelInts[state]['overlaps'].append('0')
                    iso = peakNo - c12Idx

                    # compute the mz values of the isotopic peaks according to the given charge
                    massesAndRelInts[state]['MZs'][iso] = monoMZ + (iso * neutron / charge)

                    # transfer the relative intensities to the current output isotope pattern
                    massesAndRelInts[state][relInts][iso] = self.MS1QuantIsoPatterns[state]['norm'][peakNo]

                if -1 not in massesAndRelInts[state]['MZs']:
                    massesAndRelInts[state]['overlaps'].append('0')
                    massesAndRelInts[state]['MZs'][-1] = monoMZ + (-1 * neutron / charge)
                    massesAndRelInts[state][relInts][-1] = 0

            ms1LabelStates = []

            # detect overlaps between the isotope patterns of adjacent labelling states
            # loop through the labelling states
            for labelStateA in massesAndRelInts:
                ms1LabelStates.append((labelStateA, massesAndRelInts[labelStateA]['monoMZ']))
                for labelStateB in massesAndRelInts:
                    # make sure not to compare a labelling state to itself
                    if not labelStateA == labelStateB:
                        # extract the mz values
                        patternA = sorted(massesAndRelInts[labelStateA]['MZs'].values())
                        patternB = sorted(massesAndRelInts[labelStateB]['MZs'].values())

                        # if pattern A is lighter than pattern B, but
                        if min(patternA) < min(patternB) <= max(patternA):  # doesn't work otherwise
                            for peakNoA in xrange(len(patternA)):
                                for peakNoB in xrange(len(patternB)):
                                    if abs(self.ppmDistance(patternA[peakNoA],
                                                            patternB[peakNoB])) < IsotopeData.minPeakDist:
                                        # overlap from heavier pattern (misincorporation):
                                        massesAndRelInts[labelStateA]['overlaps'][peakNoA] = '1'
                                        # overlap from lighter pattern :
                                        massesAndRelInts[labelStateB]['overlaps'][peakNoB] = '2'

            ms1LabelStates = sorted(ms1LabelStates, key=itemgetter(1))
            massesAndRelInts['byIncreasingMZ'] = []
            for i in ms1LabelStates:
                massesAndRelInts['byIncreasingMZ'].append(i[0])

            for state in self.MS1QuantIsoPatterns:
                massesAndRelInts[state]['overlaps'] = ''.join(massesAndRelInts[state]['overlaps'])
                massesAndRelInts[state]['c12index'] = self.MS1QuantIsoPatterns[state]['c12index']
                massesAndRelInts[state][sumTheo] = sum(massesAndRelInts[state][relInts].values())
            self.massesAndRelIntsByCharge[charge] = massesAndRelInts

        return self.massesAndRelIntsByCharge[charge]
