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

This file performs manipulations on the Mascot modification strings,
extracting the modifications and finding their names etc.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

import ExceptionHandler as ExHa


class MascotModifications:
    def __init__(self, modifications):
        self.MascotMods = modifications
        self.allowedOverwrites = ['Dimeth']

    def resolveDualModifiedAminoAcids(self, peptide):
        if isinstance(peptide, dict):
            variableMods = peptide['modsVariable']
            fixedMods = peptide['modsFixed']
        else:
            variableMods = peptide.modsVariable
            fixedMods = peptide.modsFixed

        combinedMods = ''
        for pos, varMod in enumerate(variableMods):
            fixMod = fixedMods[pos]
            if fixMod != '0' and varMod != '0':
                # both fixed and variable mods assigned to the same location:
                # use the variable mod by preference
                fixedMods = fixedMods[:pos] + '0' + fixedMods[pos + 1:]
                combinedMods += varMod
            elif fixMod != '0':
                combinedMods += fixMod
            else:
                combinedMods += varMod

        if isinstance(peptide, dict):
            peptide['modsFixed'] = fixedMods
            peptide['modsCombined'] = combinedMods
        else:
            peptide.modsFixed = fixedMods
            peptide.modsCombined = combinedMods
        return

    def generateModifications(self, peptide):

        self.resolveDualModifiedAminoAcids(peptide)

        sequence = peptide['peptide']
        # first deal with the variable modifications
        variableModDict = self.convertModString2Dict(peptide['modsVariable'])
        varNames = variableModDict.keys()
        varNames.sort()

        # lists for the modification data
        modStringList = []
        varModStringList = []
        peptideModsList = []
        posModDict = {}

        modStringAmbiguous = 0

        # loop through the modifications
        for name in varNames:
            number = 0
            # loop through different amino acid specificities
            aaList = []
            for amino in variableModDict[name]:
                num = len(variableModDict[name][amino])
                number += num
                # variable_modstring only needs the number, name and specificity
                varModStringList.append('%i %s (%s)' % (num, name, amino))
                for pos in variableModDict[name][amino]:
                    # positional_modstring needs number, name, specificity and location data
                    posModDict[pos] = name
                    # create mini list for the peptidemods with position and aa
                    if pos == 0:
                        aaList.append('NTerm')
                        if amino not in ['N-term', 'Protein N-term']:
                            raise ExHa.MascotModificationError('Modification %s found at N-term, specificity is %s' % (
                                name, variableModDict[name].keys()))
                    elif pos > len(sequence):
                        aaList.append('CTerm')
                        if amino != 'C-term':
                            raise ExHa.MascotModificationError('Modification %s found at C-term, specificity is %s' % (
                                name, variableModDict[name].keys()))
                    else:
                        aa = sequence[pos - 1]
                        aaList.append('%i%s' % (pos, aa))
                        if aa not in amino:
                            raise ExHa.MascotModificationError("Modification %s found at '%s', specificity is %s" % (
                                name, aa, variableModDict[name].keys()))

                # count possible sites for modifications
                if amino in ['N-term', 'C-term']:
                    numSites = 1
                else:
                    numSites = 0
                    for aa in amino:
                        numSites += sequence.count(aa)

                if numSites > num:
                    modStringAmbiguous = 1

            # peptidemods needs the number, name and list of affected AA with location
            aaList.sort()
            peptideModsList.append('%i %s:%s' % (number, name, ','.join(aaList)))

            # modstring only needs number and name
            modStringList.append('%i %s' % (number, name))

        # now analyse the fixed modifications
        fixedModDict = self.convertModString2Dict(peptide['modsFixed'])
        fixedNames = fixedModDict.keys()
        fixedNames.sort()

        fixedModStringList = []

        # loop through the modifications
        for name in fixedNames:
            number = 0
            # loop through different amino acid specificities
            for amino in fixedModDict[name]:
                num = len(fixedModDict[name][amino])
                number += num
                # fixed_modstring needs number, name and specificity
                fixedModStringList.append('%i %s (%s)' % (num, name, amino))

                for pos in fixedModDict[name][amino]:
                    # positional_modstring needs number, name, specificity and location data,
                    # accumulates with variable mod data
                    posModDict[pos] = name

        # order the positional modifications by their sequence order
        keys = posModDict.keys()
        keys.sort()
        posModStringList = []
        for pos in keys:
            name = posModDict[pos]
            posModStringList.append('%s:%i' % (name, pos))

        # turn the lists into strings and put in dictionary
        modStringsDict = dict(modstring='; '.join(modStringList), variable_modstring=';'.join(varModStringList),
                              fixed_modstring=';'.join(fixedModStringList),
                              positional_modstring=';'.join(posModStringList),
                              peptidemods='; '.join(peptideModsList), modstringambiguous=modStringAmbiguous)

        return modStringsDict

    def convertModString2Dict(self, modstring):
        # convert the string of digits into a cascading dictionary:  mod name - mod specificity -> list of locations
        locDict = {}
        for idx, mod in enumerate(modstring):
            if mod != '0':
                # modification present
                foundMod = self.MascotMods[mod]

                # create sub dictionary for modification name
                subDict = locDict.setdefault(foundMod['name'], {})
                # create sub-sub-dictionary for specificity
                subDict.setdefault(foundMod['amino'], []).append(idx)
        return locDict
