"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This holds groups of peptides together with the protein accessions of proteins
they belong to. Stats about scores and counts are based on these groups and
later also quantification

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here:
https://github.com/protcode/isob/
"""


class PeptideSet:
    def __init__(self, accession, logger):
        self.accessions = []
        self.logger = logger
        self.peptides = set()
        self.peptideData = {}
        self.ssm = 0  # Count of total peptides / spectra passing FDR filter
        self.count_novel = 0  # count of novel (will be nr /fdr ) peptides in set. Used to filter out sets when 0
        self.count_hook = 0  # count of number of hook peptides in set (will be nr / fdr).Used to filter out sets when 0
        self.score_hook = 0  # Sum of Mascot scores of all best-scoring hook peptides in set passing FDR filter
        self.score_total = 0  # Sum of Mascot scores of all best-scoring peptides in set passing FDR filter
        self.protein_group_no = 0
        self.is_reverse_hit = False
        self.reverse_protein_hits = 0
        self.upm = set()
        self.validpeptides = set()  # keeps peptide seqs which are not FDR filtered
        self.addAccession(accession)
        self.validfirstuse = 0  # flag: all peptides in set pass FDR criteria AND at least one is first use of sequence

    def __str__(self):
        accStr = ','.join(self.accessions)
        pepStr = ','.join(self.peptides)
        return '%s\n%s' % (accStr, pepStr)

    def addSequence(self, data):
        """
        @brief adds a single sequence to the peptide set
        @param data <ndarray>: containing the peptide data, includes sequence, hook score, peptide score
        @return:
        """
        peptideData = self.peptideData
        # test if the sequence has been seen before
        seq = data['sequence']
        if seq in peptideData:  # This will be the case with repeat peptides in sequences.
            # sequence already present so check the numbers
            current = peptideData[seq]
            current['number'] += int(data['numpeps'])
            # need to test each field individually
            if data['hook'] > current['isHook']:
                current['isHook'] = int(data['hook'])

            if data['hookscore'] > current['hookScore']:
                current['hookScore'] = data['hookscore']

            if data['pepscore'] > current['pepScore']:
                current['pepScore'] = data['pepscore']

            if data['bestczrank'] < current['bestCZrank']:
                current['bestCZrank'] = int(data['bestczrank'])
            peptideData[seq] = current.copy()
        else:
            # new sequence
            self.peptides.add(seq)
            peptideData[seq] = self.pepData2dict(data)

    def addAccession(self, accession):
        """
        @brief adds an accession to the list of accessions: these are the accessions with the same, full peptide set
        @param accession <string>: the accession id
        @return:
        """
        self.accessions.append(accession)
        # Mascot adds the hashes when creating the reverse or random sequences for the docoy database
        if accession.startswith('DD') or accession.startswith('###REV###') or accession.startswith('###RND###'):
            self.reverse_protein_hits += 1

        if self.reverse_protein_hits == len(self.accessions):
            self.is_reverse_hit = True
        else:
            self.is_reverse_hit = False

    @staticmethod
    def pepData2dict(pepData):
        """
        @brief converts ndarray of peptide data to a dictionary for adding to the self.peptideData dictionary
        @param pepData <ndarray>: containing the peptide data, includes sequence, hook score, peptide score
        @return:
        """
        return dict(isHook=pepData['hook'], pepScore=pepData['pepscore'], hookScore=pepData['hookscore'],
                    bestCZrank=pepData['bestczrank'], number=pepData['numpeps'], seq_start=pepData['start'],
                    seq_end=pepData['end'])

    def filterFDRpeps(self, fdrscore):
        peptideData = self.peptideData
        for seq in peptideData:
            pepobj = peptideData[seq]
            if pepobj['pepScore'] < fdrscore:
                #  No individual peptide information possible. All occurrences are below FDR threshold
                pepobj['failed_fdr_filter'] = 1
            else:
                pepobj['failed_fdr_filter'] = 0
                self.validpeptides.add(seq)

    def calcSetStatistics(self):
        """
        @brief goes through all the sequences and calculates the scores for the set
        @return:
        """
        peptideData = self.peptideData
        for seq in peptideData:
            peptide = peptideData[seq]
            if not peptide['failed_fdr_filter']:
                self.score_total += peptide['pepScore']
                self.ssm += int(peptide['number'])
                if peptide['isHook']:
                    self.score_hook += peptide['hookScore']
                    self.count_hook += 1

    def adddataFromSet(self, peptide):
        pep = self.peptideData[peptide['peptide']]
        peptide['protein_group_no'] = self.protein_group_no
        peptide['is_duplicate'] = 1
        peptide['is_first_use_of_sequence'] = pep['is_first_use_of_sequence']

        if peptide['is_unique']:
            self.upm.add(peptide['peptide'])

        if peptide['is_hook'] == pep['isHook']:
            # same hook state
            if pep['isHook'] == 1 and peptide['score'] == pep['hookScore']:
                # matches best hook score: so not duplicate
                peptide['is_duplicate'] = 0
            elif pep['isHook'] == 0 and peptide['score'] == pep['pepScore']:
                # no hook data so needs to match pepScore
                peptide['is_duplicate'] = 0
