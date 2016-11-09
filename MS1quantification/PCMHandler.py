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

class PCMdata():
    def __init__(self, mascotPeptide, MSMSheader, hasFixedLabels):
        """
        @brief creates the PTM data class with data from Mascot and pyMSsafe
        @param mascotPeptide :
        @param MSMSheader:
        @return:
        """

        self.query = mascotPeptide['query']
        self.rawsequence = mascotPeptide['sequence']
        self.sequence = '+%s-' % mascotPeptide['sequence']
        self.rank = mascotPeptide['pepno']
        self.isHook = mascotPeptide['is_hook']
        self.modsVariable = mascotPeptide['modsVariable']
        self.modsFixed = mascotPeptide['modsFixed']
        self.modsCombined = ''
        self.score = mascotPeptide['score']
        self.spec_id = MSMSheader['spec_id']
        self.mz = MSMSheader['precmz']
        self.RTmsms = MSMSheader['rt']
        self.RTapex = MSMSheader['rtapex']
        self.scanApex = MSMSheader['scanapex']
        self.charge = MSMSheader['charge']
        self.inten = MSMSheader['inten']
        self.smooth = MSMSheader['smooth']
        self.survey = MSMSheader['survey_spec']
        self.thresh = MSMSheader['thresh']
        self.PCM = '%s|%i|%s' % (self.sequence, self.charge, self.modsVariable)
        self.hasFixedLabels = hasFixedLabels

        pass

    def parseMS1QuantMods(self, allMods, ms1QuantMods):
        """
        @brief parses modstrings removing ms1 labelling modifications and finding simplest form of positionless
                  modifications
        @param varMods <dictionary>: containing the variable modification data from the Mascot search
        @param ms1QuantMods <list>: modstring identities for ms1 labelling modifications
        @return:
        """

        # remove ms1 labelling modifications from the positional modstring
        modstring = self.modsCombined

        for mod in ms1QuantMods:
            modstring = modstring.replace(mod['id'], '0')
        self.nonMS1QuantMods = modstring

        # count frequency of the modifications
        modAccumlation = {}
        modElements = {}
        for m in modstring:
            if m in allMods:
                for e in allMods[m]['elem']:
                    try:
                        modElements[e] += allMods[m]['elem'][e]
                    except KeyError:
                        modElements[e] = allMods[m]['elem'][e]

                if m in modAccumlation:
                    modAccumlation[m] += m
                else:
                    modAccumlation[m] = m

        # build the new simple mod string
        keys = modAccumlation.keys()
        keys.sort()
        simpleModString = ''
        for k in keys:
            simpleModString += modAccumlation[k]

        self.simpleModString = simpleModString
        self.modElements = modElements

        self.PCMsf = '%s|%i|%s' % (self.sequence, self.charge, self.simpleModString)
        self.PMsf = '%s|%s' % (self.sequence, self.simpleModString)

        # if self.rawsequence == 'IISNRGENSCK':
        #     pass

    def doHYPEPLEXqc(self, HYPmods, sites, labelDeltas, needed):
        """
        @brief perform differential labelling analysis for hyperplexed samples
        @param HYPmods:
        @param sites:
        @param labelDeltas:
        @return:
        """

        seq = self.sequence
        mods = self.modsVariable

        self.hyperplexQC = False
        self.ms1QuantState = None

        HYPsites = {}
        found = []
        labNTerm = set()
        labInternal = set()

        labels = {}
        for s in sites:
            labels[s] = set()

        deltaMeasured = 0

        for i, mod in enumerate(mods):
            aa = seq[i]

            if aa in sites:
                if mod in HYPmods and aa in HYPmods[mod]['amino']:
                    labels[aa].add(mod)
                    deltaMeasured += int(HYPmods[mod]['delta'] + 0.5)
                if aa in HYPsites:
                    HYPsites[aa] += 1
                else:
                    HYPsites[aa] = 1

        # need to identify mods in both events in the hyperplexing
        foundEvents = []
        for i, event in enumerate(needed):
            found = 0
            for aa in event:
                if len(labels[aa]) > 0:
                    found = 1
            foundEvents.append(found)

        if sum(foundEvents) != len(foundEvents):
            return self.hyperplexQC

        # first QC are there labels both on N-term and internally
        # for aa in labels:
        #     if len(labels[aa]) == 0 and aa in HYPsites:
        #         return self.hyperplexQC

        # first QC passed
        for state in labelDeltas:
            deltaExpect = 0
            for site in HYPsites:
                deltaExpect += HYPsites[site] * labelDeltas[state][site]

            if deltaExpect == deltaMeasured:
                self.hyperplexQC = True
                self.ms1QuantState = state
                break

        return self.hyperplexQC

    def doMS1QuantQC(self, ms1QuantMods, sites, labelDeltas):
        """
        @brief performs checks on ms1 labelling modifications
        @return:
        """
        seq = self.sequence
        mods = str(self.modsVariable)
        modsFix = str(self.modsFixed)

        if self.hasFixedLabels:
            # mix the fixed and variable mods
            for idx, fix in enumerate(modsFix):
                if fix != '0' and mods[idx] == '0':
                    mods = mods[:idx] + fix + mods[idx+1:]

        ms1LabellingSites = {}
        found = []
        labels = set()
        for i in range(len(seq)):
            aa = seq[i]
            m = mods[i]

            if aa in sites:
                # this is a site for labelling
                if aa in ms1LabellingSites:
                    ms1LabellingSites[aa] += 1
                else:
                    ms1LabellingSites[aa] = 1

                if m in ms1QuantMods:
                    found.append((aa, m, i, ms1QuantMods[m]['label'], int(ms1QuantMods[m]['delta'])))
                    labels.add(ms1QuantMods[m]['label'])

        self.ms1LabellingSites = ms1LabellingSites
        self.ms1LabellingQC = False
        self.ms1QuantState = None

        if len(labels) == 1:
            # only one ms1 labelling state
            self.ms1LabellingQC = True
            self.ms1QuantState = labels.pop()
        elif len(labels) > 1:
            # multiple ms1 labelling states: try to resolve
            if len(labelDeltas) > 2:
                # if more than two labels can be resolved into a single label
                # assuming Mascot missmatching and the label is the same for multiple locations

                deltaMeasured = 0
                for s in found:
                    deltaMeasured += s[4]

                for state in labelDeltas:
                    deltaExpect = 0
                    for site in ms1LabellingSites:
                        deltaExpect += ms1LabellingSites[site] * labelDeltas[state][site]

                    if deltaExpect == deltaMeasured:
                        self.ms1LabellingQC = True
                        self.ms1QuantState = state
                        break

        return self.ms1LabellingQC
