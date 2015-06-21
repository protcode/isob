"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This is a set of methods to manage the sets of peptides leading to protein
groups.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here:
https://github.com/protcode/isob/
"""

from peptideset import PeptideSet
import re
uniprot_gene_regexp = re.compile(' GN=(\w+) ')


class PeptideSetManager:
    def __init__(self, cfg, logger):
        self.peptideSets = {}
        self.logger = logger
        self.proteinhitdata = {}
        self.cfg = cfg
        self.pep2type = {}
        self.FDRdata = {}
        self.proteinFDRdata = {}

    def addSequenceData(self, sequenceData):
        """
        @brief iterates over the sequenceData and generates all possible accessions with all peptides
        @param sequenceData <ndarray>: containing sequence & accession data with hook scores
        @return:
        """
        peptideSets = self.peptideSets

        pep2type = {}
        for pep in sequenceData:

            peptideSets.setdefault(pep['accession'], PeptideSet(pep['accession'], self.logger)).addSequence(pep)
            if pep['accession'].startswith('DD') or pep['accession'].startswith('###REV###') or \
                    pep['accession'].startswith('###RND###'):
                accessiontype = 'REV'
            else:
                accessiontype = 'FWD'
            try:
                pep2type[pep['sequence']].add(accessiontype)
            except KeyError:
                pep2type[pep['sequence']] = {accessiontype}
        self.pep2type = pep2type

    def setFDRData(self, peptides):
        # filter rank 1s only, take only seq & score
        peptides = [(x['sequence'], x['score']) for x in peptides if x['pepno'] == 1]
        for seq, score in peptides:
            roundedscore = round(score)
            accessiontype = self.pep2type[seq]
            if accessiontype - {'REV'}:
                # ie there are FWD hits left, it's counted as a FWD hit
                # FWD hits are added to the first element of the value, REV to the second
                try:
                    self.FDRdata[roundedscore][0] += 1
                except KeyError:
                    self.FDRdata[roundedscore] = [1, 0]
            else:
                try:
                    self.FDRdata[roundedscore][1] += 1
                except KeyError:
                    self.FDRdata[roundedscore] = [0, 1]

    def setProteinGroupNo(self, protein_group_no):
        """Loops through all protein_group data and assigns protein_group_no, starting from protein_group_no"""
        my_group_no = protein_group_no
        peptideSets = self.peptideSets
        for key in peptideSets:
            pset = peptideSets[key]
            pset.protein_group_no = my_group_no
            my_group_no += 1

    def finalisePeptideSets(self, minNumHook, peptidescoreatthreshold):
        """
        @brief calculates statistics for peptidesets then creates a dictionary to group them by the number of peptides
        @param minNumHook <integer>: the minimum number of hook peptides for an acceptable protein
        @return:
        """
        peptidesets = self.peptideSets

        # remove the peptideSets with no Hook peptides, those only containing non unique peptides and those where the
        # FDR of all peptides is below the threshold
        accessions = peptidesets.keys()
        for acc in accessions:
            # remove peptides below FDR at this point.
            peptidesets[acc].filterFDRpeps(peptidescoreatthreshold)
            peptidesets[acc].calcSetStatistics()
            if peptidesets[acc].count_hook < minNumHook:
                del peptidesets[acc]
                self.logger.log.debug('will delete this peptideSet (%s): not enough hook peptides' % acc)
            elif peptidesets[acc].ssm == 0:
                self.logger.log.debug('will delete this peptideSet (%s): no valid FDR-filter-passed peptides' % acc)
                del peptidesets[acc]

        # categorise the protein data by the total number of peptides
        sets_by_length = {}
        for acc in peptidesets:
            sets_by_length.setdefault(len(peptidesets[acc].validpeptides), []).append(acc)

        self.setsByLength = sets_by_length

    def removeDuplicates(self):
        peptideSets = self.peptideSets
        setsByLength = self.setsByLength

        lengths = setsByLength.keys()
        lengths.sort(reverse=True)
        self.lengths = lengths

        for length in lengths:
            accessions = sorted(setsByLength[length])

            if len(accessions) == 1:
                # no need to filter as only one entry
                continue
            acc = accessions.pop()
            unique = [acc]
            allPeptideSets = [peptideSets[acc]]

            while accessions:
                acc = accessions.pop()
                newPeptideSet = peptideSets[acc]
                addNew = 1
                for peptideSet in allPeptideSets:
                    if newPeptideSet.validpeptides == peptideSet.validpeptides:
                        # peptide sets match so remove the sameset protein
                        del peptideSets[acc]
                        peptideSet.addAccession(acc)
                        addNew = 0
                        break
                if addNew:
                    unique.append(acc)
                    allPeptideSets.append(peptideSets[acc])

            setsByLength[length] = unique[:]

    def calcIsFirstUseOfSequence(self):

        peptideSets = self.peptideSets
        # creates the ordering required for the calculations
        ordered = []
        for acc in peptideSets:
            pepSet = peptideSets[acc]
            ordered.append((pepSet.count_hook, pepSet.score_total, len(pepSet.validpeptides), acc))

        ordered.sort(reverse=True)
        usedSequences = set()
        setNumber = 0

        for key in ordered:
            pepSet = peptideSets[key[3]]
            setNumber += 1
            pepSet.setNumber = setNumber
            validandfirstcount = 0
            self.logger.log.debug('pepsetdata %s peptides %s vs validpeptides %s'
                                  % (key[3], len(pepSet.peptideData), len(pepSet.validpeptides)))
            for seq in pepSet.peptideData:
                if seq in usedSequences:
                    pepSet.peptideData[seq]['is_first_use_of_sequence'] = 0
                else:
                    pepSet.peptideData[seq]['is_first_use_of_sequence'] = 1
                    if not pepSet.peptideData[seq]['failed_fdr_filter']:
                        validandfirstcount += 1
                    pepSet.count_novel += 1
                usedSequences.add(seq)
            self.logger.log.debug('valid and first for %s %s' % (key[3], validandfirstcount))
            pepSet.validfirstuse = validandfirstcount
        self.ordered = ordered
        self.peptide2set = usedSequences

    def calcProteinFDR(self):
        peptideSets = self.peptideSets
        # need to add the picked part of the algorithm; if we find two proteins of the same origin (REV / FWD) then pick
        # the one with the highest score and remove the other from FDR calculation
        for acc, pepset in peptideSets.iteritems():
            # get rounded top mascot score for peptides in proteinset
            max_score = max([x['pepScore'] for x in pepset.peptideData.values()])
            self.logger.log.debug('setname %s has max_score %s FP: %s' % (acc, max_score, pepset.is_reverse_hit))
            pepset.max_score = round(max_score)
            if not pepset.is_reverse_hit:
                # if our protein is a FWD (FP) hit
                try:
                    self.proteinFDRdata[max_score][0] += 1
                except KeyError:
                    self.proteinFDRdata[max_score] = [1, 0]
            else:
                try:
                    self.proteinFDRdata[max_score][1] += 1
                except KeyError:
                    self.proteinFDRdata[max_score] = [0, 1]

    def calcPeptideSetsQC(self, importData):
        stats = importData.proteinNumbers
        filtered = []
        pep2set = {}
        for key in self.ordered:
            pepSet = self.peptideSets[key[3]]
            if not pepSet.count_novel or not pepSet.validfirstuse:
                # remove the peptide set from the global list as it is common
                self.logger.log.debug('removing %s from peptideSets: no decent valid peps validfirstuse %s /'
                                      'count_novel %s' % (key[3], pepSet.validfirstuse, pepSet.count_novel))
                if not pepSet.validfirstuse:
                    self.logger.log.debug('will remove this due to no valid first use %s (cf unique)' % key[3])

                del self.peptideSets[key[3]]

                continue
            else:
                filtered.append(key)
                for pep in pepSet.peptides:
                    pep2set.setdefault(pep, set()).add(key[3])

            numHookIndex = min(5, pepSet.count_hook)
            if pepSet.is_reverse_hit:
                stats['reverse_protein_hits'][numHookIndex] += 1
            else:
                stats['forward_protein_hits'][numHookIndex] += 1

        # remove the keys from ordered
        self.ordered = filtered[:]
        self.peptide2set = pep2set

    def findUsedPeptides(self):
        """
        @brief scans through all valid sets and extracts the peptides specific to those sets. Only one peptide will
        be taken forward (from the first encountered set)
        """
        usedPeps = {}
        for key in self.ordered:
            pepSet = self.peptideSets[key[3]]
            for pep in pepSet.peptideData:

                if pep in usedPeps:
                    continue
                else:
                    usedPeps[pep] = pepSet.peptideData[pep].copy()
                    usedPeps[pep]['count_hook'] = pepSet.count_hook
                    usedPeps[pep]['score_total'] = pepSet.score_total
                    usedPeps[pep]['countPeptidesNonRed'] = len(pepSet.validpeptides)
        self.logger.log.debug('there are %s usedPeps' % len(usedPeps))
        return usedPeps

    def addPeptideSetDBdata(self, hdf5results, proteinscore2fdr):

        peptideSets = self.peptideSets
        proteinhitdata = self.proteinhitdata
        psKeys = peptideSets.keys()

        protsetIDs = range(1, len(psKeys)+1)
        proteinSetDBlist = []
        topScoringProtein = [0, None, 0]
        acc = ''
        for idx, key in enumerate(psKeys):
            psID = protsetIDs[idx]
            pepSet = peptideSets[key]
            protein_fdr = proteinscore2fdr[pepSet.max_score][3]
            self.logger.log.debug('protein_fdr %s for protein %s' % (protein_fdr, psID))
            pepSet.protein_group_no = psID
            for acc in pepSet.accessions:
                try:
                    proteindata = proteinhitdata[acc]
                    description = proteindata[0]
                    mw = proteindata[1]
                    if uniprot_gene_regexp.search(description):
                        gene_name = uniprot_gene_regexp.search(description).group(1)
                        if acc.startswith('DD') or acc.startswith('###REV###') or acc.startswith('###RND###'):
                            gene_name = '###%s###' % gene_name

                    else:
                        gene_name = None
                except KeyError:
                    self.logger.log.info('missing data for hit %s' % acc)
                    mw = 0
                    description = None
                    gene_name = None
                proteinSetDBlist.append(dict(protein_group_no=psID, protein_id=acc, mw=mw, description=description,
                                             total_score=pepSet.score_total, ssm=pepSet.ssm,
                                             hssm=pepSet.count_hook, gene_name=gene_name,
                                             upm=len(pepSet.upm),
                                             hookscore=pepSet.score_hook, is_reverse_hit=pepSet.is_reverse_hit,
                                             protein_fdr=protein_fdr, max_score=pepSet.max_score))

            #  for later we want to find the best scoring protein.
            if pepSet.score_total > topScoringProtein[0]:
                topScoringProtein = [pepSet.score_total, acc, len(pepSet.validpeptides)]

        hdf5results.writeProteinHit(proteinSetDBlist)
        return topScoringProtein

    def calcUniqueness(self, peptidescoreatthreshold=None):
        if peptidescoreatthreshold is None:
            peptidescoreatthreshold = -1
        peptideSets = self.peptideSets
        peptide2group = {}
        for a in peptideSets:

            peptides = peptideSets[a].peptides
            for peptide in peptides:

                mypeptide = peptideSets[a].peptideData[peptide]
                # only use peptide for uniqueness if its [mascot] score is greater the fdr=0.01 threshold
                if mypeptide['pepScore'] >= peptidescoreatthreshold:
                    try:
                        peptide2group[peptide].add(a)
                    except KeyError:
                        peptide2group[peptide] = {a}
        pep2unique = {}
        for p in peptide2group:

            pep2unique[p] = 1
            if len(peptide2group[p]) > 1:
                pep2unique[p] = 0
            self.logger.log.debug('for seq %s unique %s' % (p, peptide2group[p]))
        return pep2unique
