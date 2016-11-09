class ElementalCompositionCalculator:
    """
    @brief (a) calculates the elemental composition according to a peptide sequence and modifying composition
                delta_elements
           (b) applies ms1 labelling modifications to an elemental composition with respect to a provided modifying
                composition
    """

    def __init__(self, CD, ms1QuantModifications=None):
        """
        brief initialize new instance of ElementalCompositionCalculator with an optional set if fixed modifications,
                which will be applied to all elemental compositions that are computed by this instance
        @param fixedMods: <dictionary> of the form dict(C=x, H=y, O=z, ...) containing all fixed modifications
                joined together as one elemental addition
        """

        self.CD = CD
        # self.amino_acid_compositions = CD.aa_compositions_with_fixedMods
        self.amino_acid_compositions = CD.amino_acid_compositions

        self.ms1QuantModifications = ms1QuantModifications

    def calculateElementalComposition(self, sequence, variableModifications=None, noH20=False):
        """
        @brief converts the peptide sequence into an elemental composition according to the amino acid compositions in
                the class variable amino_acid_compositions
        @param sequence <string> peptide sequence to be converted
        @param variableModifications <dictionary> of the form dict(C=x, H=y, O=z, ...) containing all modifications
                joined together as one elemental addition
        @return elementalComposition <dictionary>
        """

        if noH20:
            elementalComposition = {}
        else:
            elementalComposition = {'H': 2, 'O': 1}

        # for e in self.CD.elementProperties:
        #     if not e in elementalComposition:
        #         elementalComposition[e] = 0

        for res in sequence:
            for atom in self.amino_acid_compositions[res]:
                try:
                    elementalComposition[atom] += self.amino_acid_compositions[res][atom]
                except KeyError:
                    elementalComposition[atom] = self.amino_acid_compositions[res][atom]

            #        if self.fixedMods:
            #            self.applyFixedModifications(elementalComposition)
        tmp = elementalComposition.copy()
        if variableModifications:
            self.applyVariableModsToElementalComposition(elementalComposition, variableModifications)

        return elementalComposition

    @staticmethod
    def applyVariableModsToElementalComposition(elementalComposition, modificationDefinition):
        """
        @brief applies the variable modifications stored in the instance variable self.fixedMods
        @param variableModifications <dictionary> of the form dict(C=x, H=y, O=z, ...) containing all modifications
                joined together as one elemental addition
        """
        for element in modificationDefinition:
            try:
                elementalComposition[element] += modificationDefinition[element]
            except KeyError:
                elementalComposition[element] = modificationDefinition[element]
