"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2016 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive / Fusion
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This is a basically a container for each .dat file

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""

# python modules
import math

# commonUtils modules
from CommonUtils.MascotModificationsHandler import MascotModifications

# application modules


class FileData:

    def __init__(self, filePath, hdfObj, log):

        self.filename = filePath
        self.hdfObject = hdfObj
        self.log = log
        self.score2fdr = {}
        self.proteinscore2fdr = {}
        empty = dict(mean=0, sd=0, num=0, sum=0, sumsq=0)
        self.reporterStatistics = empty.copy()
        self.precursorStatistics = dict(ppm=empty.copy(), mass=empty.copy())
        self.proteinNumbers = dict(target_protein_hits={0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0},
                                   decoy_protein_hits={0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0})

        self.acquired_spectra = 0
        self.mascot_matched_spectra = 0
        self.spectra_in_qc_proteins = 0
        self.quantified_spectra = 0
        self.numSpectraAllReporters = 0

    def addMascotModifications(self, mods):
        modDict = {}
        for m in mods:
            modDict[m['id']] = dict(name=m['name'], id=m['id'], da_delta=m['da_delta'], amino=m['amino'],
                                    relevant=m['relevant'], modtype=m['modtype'])
        self.modsHandler = MascotModifications(modDict)

    def calcStats(self, stats):
        if stats['num'] > 1:
            stats['mean'] = stats['sum'] / stats['num']
            if stats['num'] > 3:
                stats['sd'] = math.sqrt((stats['sumsq'] - stats['sum'] *
                                         stats['sum'] / stats['num']) / (stats['num'] - 1))

    def addReporterDeltas(self, deltas):

        for d in deltas:
            self.reporterStatistics['num'] += 1
            self.reporterStatistics['sum'] += d
            self.reporterStatistics['sumsq'] += d * d

    def calcReporterStats(self):
        self.calcStats(self.reporterStatistics)

    def addPrecursorDelta(self, mass, da_delta):
        ppm = da_delta / mass * 1000000

        self.precursorStatistics['ppm']['num'] += 1
        self.precursorStatistics['ppm']['sum'] += ppm
        self.precursorStatistics['ppm']['sumsq'] += ppm * ppm

        self.precursorStatistics['mass']['num'] += 1
        self.precursorStatistics['mass']['sum'] += da_delta
        self.precursorStatistics['mass']['sumsq'] += da_delta * da_delta

    def calcPrecursorStats(self):
        self.calcStats(self.precursorStatistics['ppm'])
        self.calcStats(self.precursorStatistics['mass'])

    def combineStatistics(self, newData):

        for key in ['num', 'sum', 'sumsq']:
            self.reporterStatistics[key] += newData.reporterStatistics[key]
        self.calcReporterStats()

        for tp in ['ppm', 'mass']:
            for key in ['num', 'sum', 'sumsq']:
                self.precursorStatistics[tp][key] += newData.precursorStatistics[tp][key]
        self.calcPrecursorStats()

        self.acquired_spectra += newData.acquired_spectra
        self.mascot_matched_spectra += newData.mascot_matched_spectra
        self.spectra_in_qc_proteins += newData.spectra_in_qc_proteins
        self.quantified_spectra += newData.quantified_spectra
        self.numSpectraAllReporters += newData.numSpectraAllReporters

    def getSummaryStatisticsDict(self):
        data = {}
        pn = self.proteinNumbers

        for i in range(6):
            data['%i hooks target_protein_hits' % i] = pn['target_protein_hits'][i]
            data['%i hooks decoy_protein_hits' % i] = pn['decoy_protein_hits'][i]

        data.update(self.getNumbersAccuracy())
        return data

    def getSampleDataDict(self, sample_id, fileName, runTime):
        data = dict(sample_id=sample_id,  source_file=fileName, runtime=runTime)
        data.update(self.getNumbersAccuracy())
        return data

    def getNumbersAccuracy(self):
        data = dict(mean_precursor_ion_accuracy=self.precursorStatistics['ppm']['mean'],
                    sd_precursor_ion_accuracy=self.precursorStatistics['ppm']['sd'],
                    mean_reporter_ion_accuracy=self.reporterStatistics['mean'],
                    sd_reporter_ion_accuracy=self.reporterStatistics['sd'],
                    acquired_spectra=self.acquired_spectra,
                    mascot_matched_spectra=self.mascot_matched_spectra,
                    spectra_in_qc_proteins=self.spectra_in_qc_proteins,
                    quantified_spectra=self.quantified_spectra)

        return data
