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

# CommonUtils imports
from CommonUtils.IsotopeCalculator import IsotopeCalculator

# project imports


class MS1QuantOperations:
    def __init__(self, ms1QuantModifications, CD, cfg, isoCalc=None):
        """
        @brief initialises the MS1QuantOperations class
        @param ms1QuantModifications:
        @param cfg:
        @param isoCalc:
        @return:
        """
        self.cfg = cfg
        self.ms1QuantModifications = ms1QuantModifications
        #        self.minInten = cfg.calcMinInten
        self.IC = isoCalc
        self.CD = CD
        self.amino_acid_compositions = CD.amino_acid_compositions
        self.ms1quant_amino_acid_compositions = {}
        self.ms1quant_amino_acids = ''

        self.setModifiedAmminoAcidCompositions()

    def setModifiedAmminoAcidCompositions(self):

        for sm in self.ms1QuantModifications:
            lbl = sm['label']

            if lbl not in self.ms1quant_amino_acid_compositions.keys():
                self.ms1quant_amino_acid_compositions[lbl] = {}

            modAAComp = self.ms1quant_amino_acid_compositions[lbl]

            for modifiedAA in sm['amino']:
                if modifiedAA not in self.ms1quant_amino_acids:
                    self.ms1quant_amino_acids += modifiedAA

                if modifiedAA not in modAAComp:
                    modAAComp[modifiedAA] = self.amino_acid_compositions[modifiedAA].copy()

                for element in sm['elem']:
                    try:
                        modAAComp[modifiedAA][element] += sm['elem'][element]
                    except KeyError:
                        modAAComp[modifiedAA][element] = sm['elem'][element]

    def getCoreSequence(self, sequence):
        for modifiedAA in self.ms1quant_amino_acids:
            sequence = sequence.replace(modifiedAA, '')
        return sequence

    def getComplementaryIsopattern(self, sequence, labelVersion):

        complementaryElemComp = {}
        for e in self.CD.elementProperties:
            complementaryElemComp[e] = 0

        for modifiedAA in self.ms1quant_amino_acid_compositions[labelVersion]:
            occurences = sequence.count(modifiedAA)
            modifAAComp = self.ms1quant_amino_acid_compositions[labelVersion][modifiedAA]

            if occurences:
                for element in self.ms1quant_amino_acid_compositions[labelVersion][modifiedAA]:
                    complementaryElemComp[element] += occurences * modifAAComp[element]

        if not self.IC:
            # self.IC = IsotopeCalculator(self.minInten)
            self.IC = IsotopeCalculator(self.cfg, self.CD)

        complementaryIsoPattern = self.IC.calcIsotopeDistributionFromElemComp(complementaryElemComp)

        return complementaryIsoPattern

    def applyMS1QuantModToElementalComposition(self, sequence, elementalComposition):
        """
        @brief apply the modifications defined in ms1QuantModifications to the elementalComposition and return it

        @param sequence <string> peptide sequence that is to be modified
        @param elementalComposition <dictionary> elemental composition of the peptide that is to be modified
        @return: dictionary containing an identifier for each modification and the corresponding modified elemental
                composition and names of the corresponding modifications (e.g., 'SILAC:KR+4')
        """

        if self.ms1QuantModifications is None:
            raise ValueError("No MS1 quant modifications have been specified for this instance of class "
                             "ElementalCompositionCalculator(); can be specified on initialization")

        ms1QuantModifications = self.ms1QuantModifications

        retDict = {}

        for mod in ms1QuantModifications:
            lbl = ms1QuantModifications[mod]['label']

            if lbl in retDict:
                tmpElementalComposition = retDict[lbl]['elemComposition']
                retDict[lbl]['modNames'].append(ms1QuantModifications[mod]['name'])
            else:
                tmpElementalComposition = elementalComposition.copy()
                retDict[lbl] = {'modNames': [ms1QuantModifications[mod]['name']]}

            for modifiedAA in ms1QuantModifications[mod]['amino']:
                for occurence in range(sequence.count(modifiedAA)):
                    for element in ms1QuantModifications[mod]['elem']:
                        tmpElementalComposition[element] += ms1QuantModifications[mod]['elem'][element]

            retDict[lbl]['elemComposition'] = tmpElementalComposition

        return retDict
