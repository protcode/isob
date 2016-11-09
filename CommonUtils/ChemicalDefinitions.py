import copy
import os
from pathlib import Path
from ConfigManager import ConfigManager


class ChemicalDefinitions:
    """
    @brief Container class for elemental properties and elemental compositions of amino acids.
    Fixed modifications can be applied on initialization and are then incorporated in all amino acid compositions
    accessed from this instance by the instance variable self.aa_compositions_with_fixedMods
    """
    # load monoisotopic masses, average masses, and stable isotopic abundances of all elements used from cfg
    cfgPath = Path(__file__).resolve().parent.joinpath('elements.cfg')
    cfg = ConfigManager(str(cfgPath))

    elementProperties = {}
    for e in cfg.parameters:
        elem = cfg.parameters[e]

        elemDict = dict(monoIsoMass=float(elem['monoisomass']), averageMass=float(elem['averagemass']),
                        c12idx=int(elem['c12idx']), isotopicAbundances=elem['isotopicabundances'])
        elementProperties[elem['symbol']] = elemDict.copy()

    # specify lists of elements used in standard elemental compositions and elements used for heavy isotope labelling
    labelElements = ['13C', '15N', '2H', '18O']
    standardElements = ['C', 'H', 'N', 'O', 'S', 'P', 'Se']

    #  unmodified amino acid elemental compositions
    amino_acid_compositions = {'A': {'C': 3, 'H': 5, 'N': 1, 'O': 1},  # alanine
                               'R': {'C': 6, 'H': 12, 'N': 4, 'O': 1},  # arginine
                               'N': {'C': 4, 'H': 6, 'N': 2, 'O': 2},  # asparagine
                               'D': {'C': 4, 'H': 5, 'N': 1, 'O': 3},  # aspartic acid
                               'C': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1},  # cysteine
                               'E': {'C': 5, 'H': 7, 'N': 1, 'O': 3},  # glutamic acid
                               'Q': {'C': 5, 'H': 8, 'N': 2, 'O': 2},  # glutamine
                               'G': {'C': 2, 'H': 3, 'N': 1, 'O': 1},  # glycine
                               'H': {'C': 6, 'H': 7, 'N': 3, 'O': 1},  # histidine
                               'I': {'C': 6, 'H': 11, 'N': 1, 'O': 1},  # isoluecine
                               'L': {'C': 6, 'H': 11, 'N': 1, 'O': 1},  # leucine
                               'K': {'C': 6, 'H': 12, 'N': 2, 'O': 1},  # lysine
                               'M': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1},  # methionine
                               'F': {'C': 9, 'H': 9, 'N': 1, 'O': 1},  # phenylalanine
                               'P': {'C': 5, 'H': 7, 'N': 1, 'O': 1},  # proline
                               'S': {'C': 3, 'H': 5, 'N': 1, 'O': 2},  # serine
                               'T': {'C': 4, 'H': 7, 'N': 1, 'O': 2},  # threonine
                               'W': {'C': 11, 'H': 10, 'N': 2, 'O': 1},  # tryptophan
                               'Y': {'C': 9, 'H': 9, 'N': 1, 'O': 2},  # tyrosine
                               'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1},  # valine
                               'U': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'Se': 1},  # selenocysteine
                               '+': {'C': 0, 'H': 2, 'N': 0, 'O': 0},  # N-terminus
                               '-': {'C': 0, 'H': 0, 'N': 0, 'O': 1},  # C-terminus
                               }

    averagine_composition = {'C': 4.91, 'H': 7.75, 'N': 1.38, 'O': 1.48, 'S': 0.04}

    heavy_isotopes = {'tmt': {'13C': 4, '15N': 1},
                      'itraq': {'13C': 2.25, '15N': 0.75, '18O': 0.5},
                      '8traq': {'13C': 6.5, '15N': 1.5},
                      'none': {}}

    def __init__(self, fixedMods=None):
        """
        @param fixedMods: <dictionary> of fixed modifications that are to be applied to individual amino acids
        """
        # get class variables to the local scope
        elementProperties = ChemicalDefinitions.elementProperties
        unmodified_aa_compositions = ChemicalDefinitions.amino_acid_compositions

        # modify the elemental compositions of amino acids with fixed modifications
        aa_compositions_with_fixedMods = copy.deepcopy(unmodified_aa_compositions)
        if fixedMods:
            for aa in fixedMods:
                for element in fixedMods[aa]['elem']:
                    try:
                        aa_compositions_with_fixedMods[aa][element] += fixedMods[aa]['elem'][element]
                    except KeyError:
                        aa_compositions_with_fixedMods[aa][element] = fixedMods[aa]['elem'][element]

        # provide the modified amino acids in an instance variable to make them accessible from outside the class
        self.aa_compositions_with_fixedMods = aa_compositions_with_fixedMods
