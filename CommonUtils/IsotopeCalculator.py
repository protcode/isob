from pathlib import Path
import os
import copy


class IsotopeCalculator:
    """
    @brief computes the theoretical isotope distribution for a given elemental composition
    """

    def __init__(self, cfg, CD):
        """
        @brief initialize an instance of the IsotopeCalculator class
        @param cfg <config object>: config data holding the thresholds for isotope calculation
        @param CD instance of ChemicalDefinitions, initialized with the corresponding fixed modifications
        """
        # initialise the required variables from the config data
        self.minInten = cfg.parameters['isotopecalc']['calcmininten']
        self.relIntensityFilter = cfg.parameters['isotopecalc']['reportmininten']
        self.minMissIncorporation = cfg.parameters['isotopecalc']['minmissincorp']

        # get the elemental composition of the individual amino acids, modified with the fixed modifications
        self.amino_acid_compositions = CD.amino_acid_compositions

        # get the monoisotopic mass, average mass, and isotopic abundances of the elements
        self.elementProperties = CD.elementProperties

        # get lists of the elements that unlabeled and labeled amino acids consist of
        # (for use isotope pattern construction)
        self.standardElements = CD.standardElements
        self.labelElements = CD.labelElements

        # build up a dictionary containing the isotopic patterns of each individual amino acid
        self.amino_pattern_lookup = None
        self.buildAminoPatternLookup()

        # write the dictionary the individual amino acid isotopic patterns to a file (only used for testing purposes)
        self.aminoPatternLookupFile = None

    # self.writeAminoPatternLookup()

    def buildAminoPatternLookup(self):
        """
        @brief build up a dictionary containing the isotopic patterns of each individual amino acid and store it into
                an instance variable
        """
        self.amino_pattern_lookup = {}

        for amino in self.amino_acid_compositions.keys():
            self.amino_pattern_lookup[amino] = self.calcIsotopeDistributionFromElemComp(
                self.amino_acid_compositions[amino], False)

    def writeAminoPatternLookup(self):
        """
        @brief write the dictionary the individual amino acid isotopic patterns to a file
        (only used for testing purposes)
        """
        self.aminoPatternLookupFile = Path(
            'C:\Users\holfra\Desktop\AminoPatternLookup_minInten_' + str(self.minInten).replace('.', '_') + '.tsv')

        aplOut = open(str(self.aminoPatternLookupFile), 'w')

        for i in self.amino_pattern_lookup:
            aplOut.write('%s\t' % i)

            for j in self.amino_pattern_lookup[i]['norm']:
                aplOut.write('%.20f\t' % j)

            aplOut.write('\n')
        aplOut.close()

    def filterRelIntensities(self, inPattern):
        """
        @brief filter the relative intensities in isotope pattern, keep only those that exceed a the threshold given in
                the reportMinInten variable of the MS1QuantConfig.cfg file
        @param inPattern: <dictionary> containing the isotope pattern that is to be filtered
        @return: outPattern <dictionary> containing the filtered isotope pattern
        """
        outPattern = copy.deepcopy(inPattern)

        relIntensityFilter = self.relIntensityFilter
        for normType in ['norm', 'inten']:
            outPattern[normType] = []
            if normType == 'norm':
                index = 'c12index'
            else:
                index = 'c12inten'

            for peakNo in xrange(len(inPattern[normType])):
                # only transfer the peaks that exceed relIntensityFilter threshold
                if inPattern[normType][peakNo] >= relIntensityFilter:
                    outPattern[normType].append(inPattern[normType][peakNo])

                # adjust the C12 index in the output pattern
                elif peakNo <= inPattern[index]:
                    outPattern[index] -= 1

        return outPattern

    @staticmethod
    def convoluteIsotopePatterns(isoPatternA, isoPatternB):
        """
        @brief convolute two isotope patterns
        @param isoPatternA: <list> of relative isotopic peak intensities
        @param isoPatternB: <list> of relative isotopic peak intensities
        @return: <list> of convoluted relative isotopic peak intensities
        """
        # prepare a list with the appropriate number if zeros
        convPattern = [0] * (len(isoPatternA) + len(isoPatternB))

        # do the convolution
        for i in range(len(isoPatternA)):
            for j in range(len(isoPatternB)):
                convPattern[i + j] += isoPatternA[i] * isoPatternB[j]

        return convPattern

    def addIsotopePatterns(self, isoPatternA, isoPatternB):
        """
        @brief compute the combination of two isotope patterns
        @param isoPatternA: <dictionary>
        @param isoPatternB: <dictionary>
        @return:
        """
        addPattern = {'monomass': isoPatternA['monomass'] + isoPatternB['monomass'],
                      'avgmass': isoPatternA['avgmass'] + isoPatternB['avgmass'],
                      'c12index': isoPatternA['c12index'] + isoPatternB['c12index']}

        tmpPatternA = isoPatternA['norm'][:]
        tmpPatternB = isoPatternB['norm'][:]

        addPattern['norm'] = self.convoluteIsotopePatterns(tmpPatternA, tmpPatternB)

        return self.normalizeIsotopePattern(addPattern)

    def normalizeIsotopePattern(self, inPattern):
        """
        @ brief normalization of a newly constructed isotope patterns to the C12 peak and to the most intense peak
        @param inPattern: <dictionary> the isotope pattern that is to normalized
        @return: outPattern <dictionary> the normalized isotope pattern
        """

        # set the minimum relative intensity down to which peaks are considered
        minInten = self.minInten
        minMissInc = self.minMissIncorporation

        # copy the input pattern
        outPattern = copy.deepcopy(inPattern)

        iso = inPattern['norm']
        c12Idx = inPattern['c12index']
        c12IntIdx = inPattern['c12index']

        maxint = max(iso)
        c12Int = iso[c12Idx]
        outPattern['norm'] = []
        outPattern['inten'] = []
        maxnomalsedint = 0

        for i in iso:
            normint = i / maxint
            if normint > maxnomalsedint:
                maxnomalsedint = normint

            # processed normalised data
            if normint >= minInten:
                outPattern['norm'].append(normint)
            elif maxnomalsedint < 1:
                if normint > minMissInc:
                    outPattern['norm'].append(normint)
                else:
                    c12Idx -= 1

            # process c12 based intensities
            inten = i / c12Int
            if inten >= minInten:
                outPattern['inten'].append(inten)
            elif maxnomalsedint < 1:
                if normint > minMissInc:
                    outPattern['inten'].append(inten)
                else:
                    c12IntIdx -= 1

        outPattern['c12index'] = c12Idx
        outPattern['c12inten'] = c12IntIdx

        return outPattern

    def calcIsotopeDistributionFromSequence(self, sequence):
        """
        @brief constructs an isotope pattern for a peptide from the patterns of the individual amino acids
        @param sequence: <string> peptide sequence
        @return: isoPattern <dictionary> containing the isotope pattern of the input peptide
        """

        # add H20 to the sequence to account for the peptide termini
        tmpSeq = sequence  # + 'Z'

        # set C12 index and initialize an empty list to store isotopic peak intensities in
        c12Idx = 0
        sumiso = []

        # prepare a dictionary to hold all the isotope pattern information
        isoPattern = dict(monomass=0, avgmass=0, inten=0, norm=0, c12index=None)

        # loop through the sequence
        for i in range(0, len(tmpSeq)):
            # sum up the monoisotopic and average masses of each amino acid in the sequence
            isoPattern['monomass'] += self.amino_pattern_lookup[tmpSeq[i]]['monomass']
            isoPattern['avgmass'] += self.amino_pattern_lookup[tmpSeq[i]]['avgmass']

            # nothing to convolute for the first amino acid in the sequence
            if i is 0:
                sumiso = self.amino_pattern_lookup[tmpSeq[i]]['norm']
                continue

            # for all others, cumulatively convolute their isotope patterns
            sumiso = self.convoluteIsotopePatterns(sumiso, self.amino_pattern_lookup[tmpSeq[i]]['norm'])
            c12Idx += self.amino_pattern_lookup[tmpSeq[i]]['c12index']

        # insert the convoluted data into the prepared dictionary
        isoPattern['norm'] = sumiso
        isoPattern['c12index'] = c12Idx

        # normalize the new isotope pattern
        isoPattern = self.normalizeIsotopePattern(isoPattern)

        return isoPattern

    def calcIsotopeDistributionFromElemComp(self, elemComp, clip=True):
        """
        @brief computes the theoretical isotope distribution for an elemental composition given in elemComp
        @param elemComp <dictionary> containing at least the following keys: C, H,, N, O, S, C13, N15 and the
                corresponding elemental abundances as integers
        @return: isoPattern <dictionary> of the following form:
                    monomass    -> <float> monoisotopic mass of the supplied elemental composition
                    avgmass     -> <float> average mass of the supplied elemental composition
                    inten       -> <list> of relative intensities of the elemental isotopic peaks, normalized according
                                    to the C12 peak (C12 peak = 100%)
                    norm        -> <list> of relative intensities of the elemental isotopic peaks, normalized according
                                    to the highest peak (highest peak = 100%)
                    c12index    -> <integer> index of the c12 peak
        """

        # set the minimum relative intensity down to which peaks are considered
        minInten = self.minInten

        # prepare a dictionary to hold all the isotope pattern information
        isoPattern = dict(monomass=0, avgmass=0, inten=0, norm=0, c12index=None)

        # initialize an empty list to store isotopic peak intensities in
        sumiso = 0
        # add the labeled isotope data
        c12Idx = 0

        stdElemComp = [e for e in elemComp if e not in self.labelElements and elemComp[e] > 0]
        labelElemComp = [e for e in elemComp if e in self.labelElements and elemComp[e] > 0]

        # loop through the light isotopes of all possible elements
        for element in stdElemComp:
            # construct the isotopic distribution of the current element
            elemIsoPattern = self.elementProperties[element]['isotopicAbundances'][:]
            # elemIsoPattern.insert(0, 1)

            # loop through the abundance of the current element in the given composition
            for j in xrange(elemComp[element]):
                # cumulate the elements' monoisotopic and average masses
                isoPattern['monomass'] += self.elementProperties[element]['monoIsoMass']
                isoPattern['avgmass'] += self.elementProperties[element]['averageMass']
                c12Idx += self.elementProperties[element]['c12idx']

                # nothing to convolute for the first atom in the composition
                if sumiso == 0:
                    sumiso = elemIsoPattern
                    continue

                # for all others, cumulatively convolute their isotope patterns
                sumiso = self.convoluteIsotopePatterns(sumiso, elemIsoPattern)

                # for speed-up: clip off all isotopic peaks that are smaller than the minimum relative intensity
                # down to which peaks should be considered
                if clip:
                    for i in range(len(sumiso) - 1, 0, -1):
                        if sumiso[i] > minInten:
                            break
                    sumiso = sumiso[:i + 1]

        if sumiso == 0:
            sumiso = [1]

        # loop through the heavy isotopes used for labelling
        for element in labelElemComp:
            # loop through the abundance of the current heavy isotope in the given composition
            for j in range(elemComp[element]):
                # add the isotope's monoisotopic and average masses
                isoPattern['monomass'] += self.elementProperties[element]['monoIsoMass']
                isoPattern['avgmass'] += self.elementProperties[element]['averageMass']

                # add a leading '0' to the isotope pattern that has been constructed so far and adjust the C12 iindex
                sumiso = [0] + sumiso
                c12Idx += 1

                # convolute the heavy isotope's pattern with the pattern constructed so far
                for i in range(len(sumiso) - 1):
                    sumiso[i] += sumiso[i + 1] * self.elementProperties[element]['isotopicAbundances'][0]

        # insert the convoluted data into the prepared dictionary
        isoPattern['norm'] = sumiso[:]
        isoPattern['c12index'] = c12Idx
        isoPattern['c12inten'] = c12Idx

        # normalize the constructed isotope pattern
        isoPattern = self.normalizeIsotopePattern(isoPattern)

        return isoPattern

    def calcMassesFromElemComp(self, elemComp):
        masses = dict(monomass=0, avgmass=0, elements=elemComp.copy())

        for element in elemComp:
            masses['monomass'] += self.elementProperties[element]['monoIsoMass'] * elemComp[element]
            masses['avgmass'] += self.elementProperties[element]['averageMass'] * elemComp[element]

        return masses
