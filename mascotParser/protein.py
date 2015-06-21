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

This file is the protein container for Mascot protein information.
This holds the peptide data associated with the protein.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""


class proteinSet():
    def __init__(self, accession, peptide, pepData):
        self.accession = accession
        self.allPeptides = {peptide}
        self.peptideData = {peptide: pepData.copy()}
        self.numPeptides = 1
        if pepData['is_hook'] == 1:
            self.hookPeptides = {peptide}
            self.numHookPeptides = 1
        else:
            self.hookPeptides = set()
            self.numHookPeptides = 0

    def addPeptide(self, peptide, pepData):
        allPeps = self.allPeptides

        if peptide in allPeps:
            oldPepData = self.peptideData[peptide]

            if oldPepData['is_hook'] < pepData['is_hook']:
                oldPepData = pepData.copy()
                self.hookPeptides.add(peptide)
                self.numHookPeptides += 1
                self.numPeptides += 1
            elif oldPepData['is_hook'] == pepData['is_hook'] and oldPepData['max_score'] < pepData['max_score']:
                oldPepData = pepData.copy()
                self.numPeptides += 1

        else:
            self.peptideData[peptide] = pepData.copy()
            allPeps.add(peptide)
            self.numPeptides += 1
            if pepData['is_hook'] == 1:
                self.hookPeptides.add(peptide)
                self.numHookPeptides += 1
