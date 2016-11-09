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
from CommonUtils.ElementalCompositionCalculator import ElementalCompositionCalculator
from CommonUtils.ChemicalDefinitions import ChemicalDefinitions
from CommonUtils.averagine import AveragineModel

# project imports
from IsotopeData import IsotopeData
from MS1QuantOperations import MS1QuantOperations


class IsotopeCollection():
    """
    @brief central storage and lookup class for mapping PMs to their corresponding elemental compositions modified for
    MS1 quanitification and isotope patterns which are stored in IsotopeData() instances
    """

    def __init__(self, cfg, hyperplexed, fixedMods=None, ms1QuantModifications=None, minInten=None, isoCalcMethod='bySEQ'):
        if not isoCalcMethod in ['byEC_noClip', 'byEC_clip', 'bySEQ']:
            raise ValueError("Parameter isoCalcMethod must be either 'byEC_noClip', 'byEC_clip', or 'bySEQ'")

        self.cfg = cfg
        self.hyperplexed = hyperplexed
        self.isotopeLookup = {}
        self.CD = ChemicalDefinitions(fixedMods)
        self.ECC = ElementalCompositionCalculator(self.CD, ms1QuantModifications)
        self.IC = IsotopeCalculator(cfg, self.CD)
        self.SO = MS1QuantOperations(ms1QuantModifications, self.CD, cfg, self.IC)
        self.countPMs = 0
        self.uniquePMs = 0
        self.ms1QuantModifications = ms1QuantModifications
        self.fixedMods = fixedMods
        self.isoCalcMethod = isoCalcMethod
        self.averagineModel = AveragineModel(cfg.parameters['runtime']['datadir'], 300, 0.01, 'none', 2, cfg)

        if self.cfg.parameters['models']['primary_model'].lower() == 'exact':
            self.primary = 'exact'
            self.exact = 'primary'
            self.secondary = 'averagine'
        else:
            self.primary = 'averagine'
            self.secondary = 'exact'
            self.exact = 'secondary'

        if not self.cfg.parameters['models']['compare']:
            self.secondary = 0


    def addLookupEntry(self, pcm):
        """
        @brief Adds an entry into the isotopeLookup dictionary. They key of this entry is the PMsf (ms1-quant-free
            Peptide + modifications) string, which consists of the peptide sequence and a modification string which
            does not contain ms1-quant modifications. The value of the dictionary entry is an instance of the
            IsotopeData class, which stores elemental compositions and isotope patterns for an elemental composition
            and its ms1-quant modifications. If the isotopeLookup dictionary already contains an entry with the
            considered PMsf as its key, nothing is done.
        @param pcm: <PCMdata> containing the considered PCM information
        """
        self.countPMs += 1

        if pcm.PMsf in self.isotopeLookup.keys():
            return

        self.uniquePMs += 1
        elemComp = self.ECC.calculateElementalComposition(pcm.sequence, pcm.modElements, True)
        tmpIsoData = IsotopeData(elemComp, pcm.sequence, self.ms1QuantModifications, self.ECC, self.IC, self.SO,
                                 self.exact, self.isoCalcMethod, varMods=pcm.modElements)

        self.isotopeLookup[pcm.PMsf] = tmpIsoData


    def getMS1quantModifiedIsoPatterns(self, pcm):
        """
        @brief Returns the mz values and corresponding relative intensities for the isotopic patterns of all labelling
            states of a PM.
        @param pcm: <PCMdata> containing the considered PCM information
        @return:
        """
        if not pcm.PMsf in self.isotopeLookup:
            raise KeyError('The PM ' + pcm.PMsf + ' is not contained in this IsotopeCollection\'s isotopeLookup')

        pcmData = self.isotopeLookup[pcm.PMsf]
        # this sets the exact model data with either 2 or 1 as model name

        pcmData.getMassesAndRelInts(pcm.charge)

        for state in pcmData.ms1QuantStates:
            label = pcmData.massesAndRelIntsByCharge[pcm.charge][state]
            if self.cfg.parameters['models']['primary_model'] == 'averagine' or self.cfg.parameters['models']['compare']:
                averagineIso = self.averagineModel.findIsotopeSet(pcmData.MS1QuantIsoPatterns[state]['monomass'])
            else:
                averagineIso = dict(RelInts=dict(), inten=list(), sumTheo=0)
            label['averagine'] = averagineIso
            if self.exact == 'primary':
                label['primaryModel'] = 'exact'
                label['secondaryModel'] = 'averagine'
                label['secondary'] = averagineIso['RelInts']
                label['secondary_sumTheo'] = averagineIso['sumTheo']
            else:
                label['primaryModel'] = 'averagine'
                label['secondaryModel'] = 'exact'
                label['primary'] = averagineIso['RelInts']
                label['primary_sumTheo'] = averagineIso['sumTheo']

        return pcmData.massesAndRelIntsByCharge[pcm.charge]
