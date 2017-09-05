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

This file is the data container for Mascot peptide data, controlling the
extraction of peptide data and creation of modification data.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""


import re
pepbits = [('misscleave', 'i'),
           ('mass', 'f'),
           ('da_delta', 'f'),
           ('numionsmatched', 'i'),
           ('sequence', 's'),
           ('peaks1', 'i'),
           ('modsVariable', 's'),
           ('score', 'f'),
           ('seriesfound', 's'),
           ('peaks2', 'i'),
           ('peaks3', 'i')]

re_reverse = re.compile('(###R[A-Z]{2}###)')
re_prbits = re.compile('^(.+):\d+:(\d+):(\d+):\d+$')

class Peptide:
    def __init__(self, inputdic, datfile):
        """
        @brief initialises peptide object by parsing the data in inputdic
        @param inputdic <dictionary>: containing the mascot peptide data
        @param datfile <datfile object>: datfile parsed data and methods
        """

        # application object data
        self.cfg = datfile.cfg
        self.logs = datfile.logs
        self.datfileobj = datfile
        self.is_hook = 0
        self.useinprot = 1
        self.retain = 0
        self.czrank = 0

        for key in inputdic:
            keybits = key.split('_', 2)
            if len(keybits) == 2:
                # main part of the peptide
                self.query = int(keybits[0][1:])
                self.pepno = int(keybits[1][1:])
                self.parsePeptideString(inputdic[key])
            elif keybits[2] == 'terms':
                # catch the flanking amino acids
                self.flanks = inputdic[key].split(':')
            elif keybits[2] == 'et_mods':
                # catch the error tolerant modifications
                fields = inputdic[key].split(',')
                self.etMods = dict(da_delta=float(fields[0]), nl=float(fields[1]), desc=fields[2], name=fields[0])
                x = 1
            elif keybits[2] == 'subst':
                # amino acid substitutions
                self.substitutions = inputdic[key].split(':')
            elif keybits[2] == 'primary_nl':
                # primary neutral losses
                self.primary_nl = inputdic[key]
            else:
                # catch other data lines here and report
                self.logs.datlog.debug('Skipped unexpected data line: %s=%s' % (key, inputdic[key]))

    def doSplitProteinIDs(self, prostring):
        # essentially splits proteinlist on the ',' string
        # but because some definitions also contain commas we have to
        # re-concatenate the separate strings
        g = prostring.split(',')
        y = 0
        newlist = list()
        while y < len(g):
            testline = g[y]
            if not testline:
                break
            while testline.count(':') < 4:
                y += 1
                testline += ',' + g[y]
                testline.count(':')
            newlist.append(testline)
            y += 1
        return newlist

    def parsePeptideString(self, pepstr):
        """
        @brief parses the main peptide string
        @param pepstr <string>: data string from mascot peptide
        """
        # first split the peptide and protein sections, and split the parts

        pepprot = pepstr.split(';"')
        pepdata = pepprot[0].split(',')
        protdata = self.doSplitProteinIDs('"'+''.join(pepprot[1:]))

        # now process the peptide section
        for i, val in enumerate(pepdata):
            if pepbits[i][1] == 'i':
                # integer data
                self.__dict__[pepbits[i][0]] = int(val)
            elif pepbits[i][1] == 'f':
                # float data
                self.__dict__[pepbits[i][0]] = float(val)
            elif pepbits[i][1] == 's':
                # string data
                self.__dict__[pepbits[i][0]] = val

        prots = []
        for pr in protdata:
            if re_prbits.match(pr):
                start = re_prbits.match(pr).group(2)
                stop = re_prbits.match(pr).group(3)
                full_accession = re_prbits.match(pr).group(1).replace('"','')
                items = full_accession.split('|')
            else:
                prbits = pr.split(':')
                items = prbits[0][1:-1].split('|')
                start = prbits[2]
                stop = prbits[3]
            if len(items) > 2:
                accession = items[1]
            else:
                accession = items[0]
            if re_reverse.search(items[0]):
                rev_tag = re_reverse.search(items[0]).group(1)
                accession = rev_tag + accession

            prots.append(dict(accession=accession, start=int(start),  end=int(stop)))

        self.proteins = prots
        self.modsVariable = self.modsVariable.lower()

        self.modsFixed, self.modsRelevant = self.calcModStrings()

    def isValidSequence(self):
        """
        @brief tests the sequence against valid amino acids
        """

        invalaa = [aa for aa in self.sequence if aa not in self.cfg.parameters['general']['allowedamino']]
        return not invalaa

    def calcModStrings(self):
        """
        @brief calculates the modstings from the sequences etc
        @param mods <dictionary>: containing the modification data
        """
        mods = self.datfileobj.mods
        modsVariable = self.modsVariable
        seq = self.sequence
        pmax = len(seq)

        variable = {}
        fixed = {}
        fixedModString = ''
        relevantModString = ''

        # extract all mods
        for p, m in enumerate(modsVariable):
            # calculat the relevant modstring
            if m == 'X':
                # this is an ET modification
                relevantModString += m
                variable[p] = self.etMods
            elif m != '0':
                # variable modification found
                variable[p] = mods[m]
                if mods[m]['relevant']:
                    relevantModString += m
                else:
                    relevantModString += '0'
            else:
                relevantModString += '0'

            # calculate the fixed modstring
            if p == 0:
                # N-term: check for mods
                if '+' in mods:
                    fixedModString += '+'
                else:
                    fixedModString += '0'
            elif p > pmax:
                # C-term: check for mods
                if '-' in mods:
                    fixedModString += '-'
                else:
                    fixedModString += '0'
            elif seq[p - 1] in mods:
                # fixed mod found
                fixed[p - 1] = mods[seq[p - 1]]
                fixedModString += seq[p - 1]
            else:
                fixedModString += '0'
        try:
            i = int(relevantModString)
        except:
            pass
        return fixedModString, relevantModString


class ETpeptide(Peptide):
    def __init__(self, inputdic, datfile):
        """
        @brief initialises peptide object by parsing the data in inputdic
        @param inputdic <dictionary>: containing the mascot peptide data
        @param datfile <datfile object>: datfile parsed data and methods
        """

        # application object data
        self.cfg = datfile.cfg
        self.datfileobj = datfile
        self.is_hook = 0
        self.useinprot = 1
        self.retain = 0
        self.czrank = 0

        # find the line identifiers for the peptide
        keys = inputdic.keys()
        keys.sort()
        qryPep = keys.pop(0)

        # extract the query and peptide
        bits = qryPep.split('_')
        self.query = int(bits[0][1:])
        self.pepno = int(bits[1][1:])

        # check if there is flanking amino acid data
        terms = qryPep + '_terms'
        if terms in keys:
            keys.remove(terms)
            self.flanks = inputdic[terms].split(':')

        # check for error tolerant modifications
        etmods = qryPep + '_et_mods'
        if etmods in keys:
            keys.remove(etmods)
            fields = inputdic[etmods].split(',')
            self.etMods = dict(da_delta=float(fields[0]), nl=float(fields[1]), desc=fields[2], name=fields[0])

        # parse the main peptide data
        self.parsePeptideString(inputdic[qryPep])

        for key in keys:
            # catch other data lines here and report
            self.logs.datlog.debug('Skipped unexpected data line: %s=%s' % (key, inputdic[key]))

    def calcModStrings(self):
        """
        @brief calculates the modstings from the sequences etc
        @param mods <dictionary>: containing the modification data
        """
        mods = self.datfileobj.mods
        modsVariable = self.modsVariable
        seq = self.sequence
        pmax = len(seq)

        variable = {}
        fixed = {}
        fixedModString = ''
        relModString = ''

        # extract all mods
        for p, m in enumerate(modsVariable):
            # calculate the relevant modstring
            if m == 'X':
                relModString += m
                variable[p] = 'X'
            elif m != '0':
                # variable modification found
                variable[p] = mods[m]
                if mods[m]['relevant']:
                    relModString += m
                else:
                    relModString += '0'
            else:
                relModString += '0'

            # calculate the fixed modstring
            if p == 0:
                # N-term: check for mods
                if '+' in mods:
                    fixedModString += '+'
                else:
                    fixedModString += '0'
            elif p > pmax:
                # C-term: check for mods
                if '-' in mods:
                    fixedModString += '-'
                else:
                    fixedModString += '0'
            elif seq[p - 1] in mods:
                # fixed mod found
                fixed[p - 1] = mods[seq[p - 1]]
                fixedModString += seq[p - 1]
            else:
                fixedModString += '0'
        try:
            i = int(relModString)
        except:
            pass
        return fixedModString, relModString
