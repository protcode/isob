"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This is an API to the MS.hdf5 files and contains methods to read
data from it

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/cellzome/isobarquant
"""
# python modules

import sys
sys.path.insert(0, '..')
# commonUtils modules
from CommonUtils.hdf5Base import hdf5Base
from CommonUtils.QuantMethodHandler import QuantMethods


class HDFfile:
    def __init__(self, hdfFilePath, logs):
        self.hdfFilePath = hdfFilePath
        if hasattr(logs, 'hdf5log'):
            self.log = logs.hdf5log
        else:
            logs.setLogName('hdf5log', 'hdf5file')
            self.log = logs.hdf5log

        hdf = hdf5Base(hdfFilePath)
        hdf.readOpen()
        # get all imported data from HDF5 file
        allimports = hdf.readTable('/imports')
        # sort by dat when search run
        imports = sorted(allimports, key=lambda x: x[1], reverse=True)
        # use the lastest files
        self.importPath = imports[0]['name']
        # find the quant method
        isotopes = hdf.readTable('/rawdata/isotopes')
        if len(isotopes) > 0:
            quantMethodID = isotopes[0]['method_id']
            qm = QuantMethods()
            self.quantMethod = qm.getMethodByID(quantMethodID)
            self.quantLoc = '/rawdata/quan'

        else:
            self.quantMethod = None
            self.quantLoc = None
            self.quantExtraLoc = None
        # expect quant tables to be in rawdata - except MS1 may move to mascot area

        hdf.close()
        self.hdf = hdf

    def replaceTable(self, table, data, definition):
        self.hdf.appendOpen()
        self.hdf.replaceTable(self.importPath, table, data, definition)
        self.hdf.close()

    def readSeq2accTable(self):

        self.hdf.readOpen()
        s2a = self.hdf.readTable('/%s/seq2acc' % self.importPath)
        self.hdf.close()

        return s2a

    def readProteinTable(self):
        self.hdf.readOpen()
        proteins = self.hdf.getTable('/%s/proteins' % self.importPath)

        accession2proteininfo = {}
        for row in proteins:

            accession = row[0]
            mw = row[1]
            description = row[2]
            accession2proteininfo[accession] = (description, mw)
        self.hdf.close()
        return accession2proteininfo

    def readMascotConfigTable(self):

        self.hdf.readOpen()
        config = self.hdf.readTable('/%s/config' % self.importPath)
        self.hdf.close()

        return config

    def readModsTable(self):

        self.hdf.readOpen()
        mods = self.hdf.readTable('/%s/mods' % self.importPath)
        self.hdf.close()

        return mods

    def readStatistics(self):

        self.hdf.readOpen()
        stats = self.hdf.readTable('/%s/statistics' % self.importPath)

        statsDict = {}
        for s in stats:
            if s['statistic'][:3] == 'num':
                statsDict[s['statistic']] = int(s['value'])
            else:
                statsDict[s['statistic']] = s['value']
        self.hdf.close()

        return statsDict

    def readSpectraTable(self):
        self.hdf.readOpen()
        spectra = self.hdf.readTable('/%s/peptides' % self.importPath)
        self.hdf.close()
        return spectra

    def readMSMSheader(self):
        self.hdf.readOpen()
        data = self.hdf.readTable('/rawdata/msmsheader')
        self.hdf.close()

        return data

    def readQuan(self):
        self.hdf.readOpen()
        data = self.hdf.readTable(self.quantLoc)
        self.hdf.close()

        return data

    def readQuanExtra(self):
        self.hdf.readOpen()
        data = self.hdf.readTable(self.quantExtraLoc)
        self.hdf.close()

        return data

    def getAcquisitionTime(self):
        self.hdf.readOpen()
        acqTime = self.hdf.getDataEqual('/rawdata/parameters', 'parameter', 'MS Run Time (min)')
        self.hdf.close()
        return float(acqTime[0]['value'])

    def getTimeAndActivation(self):
        idAct = None
        quanAct = None
        self.hdf.readOpen()
        acqTime = self.hdf.getDataEqual('/rawdata/parameters', 'parameter', 'MS Run Time (min)')
        activation = self.hdf.getDataEqual('/rawdata/parameters', 'parameter', 'Activation Type')
        if len(activation) == 0:
            activation = self.hdf.getDataEqual('/rawdata/parameters', 'parameter', 'activation')
        self.hdf.close()

        actTypes = set()
        for a in activation.flat:
            actTypes.add(a['value'].upper())

        if len(actTypes) == 1:
            idAct = actTypes.pop()
            quanAct = idAct
        elif len(actTypes) == 2:
            act = actTypes.pop()
            if act == 'HCD':
                quanAct = 'HCD'
                idAct = actTypes.pop()
            else:
                quanAct = 'HCD'
                idAct = act
        else:
            pass
        return float(acqTime[0]['value']), idAct, quanAct

    def readPeptides(self):
        self.hdf.readOpen()
        peps = self.hdf.readTable('/%s/peptides' % self.importPath)
        self.hdf.close()
        return peps

    def getRelevantPeptides(self, usedPeps, hdfFileData):
        self.hdf.readOpen()

        peps = self.hdf.getDataEqual('/%s/peptides' % self.importPath, 'useinprot', 1)
        self.hdf.close()

        # filter the peptide data to only those peptides in QCed proteins
        pepDict = {}
        for pep in peps:
            if pep['sequence'] in usedPeps:
                # pep is used
                pepData = self.extractPepDataFromUsedPeps(pep)
                pepData['failed_fdr_filter'] = usedPeps[pep['sequence']]['failed_fdr_filter']
                pepData['seq_start'] = usedPeps[pep['sequence']]['seq_start']
                pepData['seq_end'] = usedPeps[pep['sequence']]['seq_end']
                hdfFileData.addPrecursorDelta(pep['mass'], pep['da_delta'])

                pepDict.setdefault(pep['query'], []).append(pepData.copy())

        hdfFileData.calcPrecursorStats()
        return pepDict

    def readQueries(self):
        self.hdf.readOpen()
        queries = self.hdf.readTable('/%s/queries' % self.importPath)
        self.hdf.close()
        return queries

    def readQueryDict(self, keys):

        queries = self.readQueries()

        queryDict = {}
        for q in queries:
            tmp = {}
            for key in keys:
                tmp[key] = q[key]

            queryDict[q['query']] = tmp.copy()

        return queryDict

    def getRelevantQueries(self, pepDict):
        queries = self.readQueries()
        queryDict = {}
        for q in queries:
            if q['query'] in pepDict:
                queryDict[q['query']] = dict(query=q['query'], spec_id=q['spec_id'],
                                             msms_id=q['msms_id'],
                                             neutral_mass=q['prec_neutmass'], delta_seq=q['delta_seq'],
                                             delta_mod=q['delta_mod'])
        return queryDict

    @staticmethod
    def extractPepDataFromUsedPeps(pep):
        return dict(query=int(pep['query']), rank=int(pep['pepno']),
                    peptide=pep['sequence'], is_hook=int(pep['is_hook']),
                    modsVariable=pep['modsVariable'], modsFixed=pep['modsFixed'],
                    mw=pep['mass'], da_delta=pep['da_delta'], score=pep['score'],
                    missed_cleavage_sites=int(pep['misscleave']))

    def getNumbers(self):
        self.hdf.readOpen()
        data = self.hdf.getDataEqual('/%s/statistics' % self.importPath, 'statistic', 'numspectra')
        acquired_spectra = int(data['value'])

        data = self.hdf.getDataEqual('/%s/statistics' % self.importPath, 'statistic', 'numspectra_nopeps')
        numNoPeps = int(data['value'])

        table = self.hdf.tableTOC['/rawdata/isotopes']
        numIsotopes = table.nrows

        data = self.hdf.getDataEqual('/rawdata/parameters', 'parameter', 'MS Run Time (min)')
        runTime = float(data[0]['value'])

        self.hdf.close()

        return acquired_spectra, acquired_spectra - numNoPeps, numIsotopes, runTime

    def readFragDeltas(self):
        self.hdf.readOpen()
        ppmDeltaData = []
        daDeltaData = []
        tablePath = '/%s/fragdelta' % self.importPath

        if tablePath in self.hdf.tableTOC:
            ppmData = self.hdf.getDataEqual('/%s/fragdelta' % self.importPath, 'errortype', 'ppm')
            daData = self.hdf.getDataEqual('/%s/fragdelta' % self.importPath, 'errortype', 'da')
        else:
            ppmData = []
            self.log.info('fragdelta table not found in %s' % self.importPath)
            daData = []

        for ppm in ppmData:
            ppmDeltaData.append([ppm['score'], ppm['idtype'], ppm['error']])
        for da in daData:
            daDeltaData.append([da['score'], da['idtype'], da['error']])
        self.hdf.close()
        return ppmDeltaData, daDeltaData

    def readImporterData(self, usedPeps, hdf):
        # read the relevant query data
        self.log.debug('reading HDF5 peptide data')
        peptides = self.getRelevantPeptides(usedPeps, hdf)
        self.log.debug('reading HDF5 query data')
        queryDict = self.getRelevantQueries(peptides)
        self.log.debug('reading HDF5 MS/MS headers data')
        headerArray = self.readMSMSheader()
        self.log.debug('reading HDF5 quantification data')

        if self.quantMethod:
            quanArray = self.readQuan()
        else:
            quanArray = None

        return peptides, queryDict, headerArray, quanArray
