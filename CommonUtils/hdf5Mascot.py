# -*- coding: utf-8 -*-
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

This file handles the hdf5 interface for Mascot derived data.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob

"""
from hdf5Base import hdf5Base
import datetime
import unicodedata

import ExceptionHandler as ExHa


class HDF5Mascot:
    def __init__(self, hdfFilePath=None, hdfBaseObj=None):

        if not hdfFilePath and not hdfBaseObj:
            raise (ExHa.HDF5consistancyError, "Neither hdfFilePath nor hdfBaseObj provided.")

        elif hdfFilePath and hdfBaseObj:
            if hdfBaseObj.filePath == hdfFilePath:
                self.hdf = hdfBaseObj
                self.hdfFilePath = hdfFilePath

            else:
                raise (ExHa.HDF5consistancyError,
                       "The provided hdfFilePath and hdfBaseObj do not reference the same file.")

        elif hdfBaseObj and not hdfFilePath:
            self.hdf = hdfBaseObj
            self.hdfFilePath = hdfBaseObj.filePath

        elif hdfFilePath and not hdfBaseObj:
            self.hdfFilePath = hdfFilePath
            self.hdf = hdf5Base(hdfFilePath)

        self.numberImports = 0
        self.importGroup = ''
        self.impRow = -1

        self.datFileName = None

    def appendOpen(self):
        self.hdf.appendOpen()

    def readOpen(self):
        self.hdf.readOpen()

    def close(self):
        self.hdf.close()

    def checkDatFilePresence(self, datFileName):
        """
        @brief checks the named dat file has been imported
        @param datFileName <string>: the datfile name including extension
        @return importGroup <string>: the name of the group that the dat file is imported to, if any
        """
        self.datFileName = datFileName

        # if '/imports' in self.hdf.groupTOC:
        if '/imports' in self.hdf.tableTOC:
            imports = self.hdf.readTable('/imports')

            self.numberImports = len(imports)

            for row, imp in enumerate(imports):
                if imp['datfile'] == datFileName:
                    self.importGroup = imp['name']
                    self.impRow = row
                    break

            return self.importGroup
        else:
            return False

    def getLatestDatfile(self):

        if '/imports' in self.hdf.tableTOC:
            imports = self.hdf.readTable('/imports')

            self.numberImports = len(imports)
            sorted(imports, key=lambda x: x['date'])
            for row, imp in enumerate(imports):
                self.importGroup = imp['name']
                self.impRow = row

            return self.importGroup
        else:
            return False

    def getPeptides(self):
        hdf = self.hdf
        tablePath = '/%s/%s' % (self.importGroup, 'peptides')
        peptides = hdf.readTable(tablePath)
        return peptides

    def getPeptidefromMSMSID(self, msms_id, rank=None):
        hdf = self.hdf
        tablePath = '/%s/%s' % (self.importGroup, 'queries')
        queryData = hdf.getDataGeneral(tablePath, "spec_id==%s" % msms_id)[0]
        query = queryData['query']

        tablePath = '/%s/%s' % (self.importGroup, 'peptides')

        peptideData = hdf.getDataGeneral(tablePath, "query==%s" % query)
        if rank:
            peptideData = [x for x in peptideData if x[1] == rank]

        if peptideData:
            return peptideData[0]
        else:
            return []

    def getQueryno2msmsid(self):
        hdf = self.hdf
        tablePath = '/%s/%s' % (self.importGroup, 'queries')
        returndict = dict()
        for x in hdf.readTable(tablePath):
            query = x['query']
            spec_id = x['spec_id']
            returndict[query] = spec_id
        return returndict

    def getMsmsid2Queryno(self):
        hdf = self.hdf
        tablePath = '/%s/%s' % (self.importGroup, 'queries')
        returndict = dict()
        for x in hdf.readTable(tablePath):
            query = x['query']
            spec_id = x['spec_id']
            returndict[spec_id] = query
        return returndict

    def deleteMascotImport(self, importGroup):
        """
        @brief delete all Mascot data from an import
        @param importGroup <string>: defining the import group name
        @return:
        """
        hdf = self.hdf

        # remove tables
        tables2delete = [x for x in hdf.tableTOC if importGroup in x]
        for table in tables2delete:
            hdf.removeTable(table)

        # remove group
        hdf.removeGroup(importGroup)

        # remove imports table entry
        if self.numberImports > 1:
            hdf.removeRows('/imports', {'name': importGroup})
        else:
            hdf.removeTable('/imports')
        self.numberImports -= 1

        return

    def deleteAllMascotImports(self, isPreCal):
        """
        @brief deletes all mascot imports
        @return:
        """
        if '/imports' in self.hdf.tableTOC:
            imports = self.hdf.readTable('/imports')

            for imp in imports:
                if isPreCal == 0 and imp['name'] == 'calibration':
                    continue

                self.deleteMascotImport(imp['name'])
            if isPreCal == 1:
                self.hdf.removeTable('/imports')
        return

    def createTables(self, importGroup, isPreCal):
        """
        @brief creates all the required tables for the importing of Mascot data
        @param importGroup <string>: the group name to be created to hold the data
        @return:
        """

        hdf = self.hdf
        # add entry to the imports table
        if '/imports' not in hdf.tableTOC:
            hdf.createTable('/', 'imports', 'Imports')
        t = datetime.datetime.now()
        strTime = t.strftime('%Y-%m-%d %H:%M:%S')
        self.impRow = len(self.hdf.tableTOC['/imports'])
        hdf.appendRows('/imports', [dict(name=importGroup, date=strTime, datfile=self.datFileName, status='running')])

        hdf.createGroup(importGroup)
        if self.importGroup != importGroup:
            self.importGroup = importGroup
        if isPreCal:
            tables = [('parameters', 'Parameter'), ('masses', 'Mass'), ('mods', 'Mod'), ('calmasses', 'CalMasses')]
        else:
            tables = [('parameters', 'Parameter'), ('masses', 'Mass'), ('mods', 'Mod'), ('queries', 'Query'),
                      ('peptides', 'Peptide'), ('proteins', 'Protein'), ('seq2acc', 'Seq2Acc'), ('index', 'Index'),
                      ('statistics', 'Statistics'), ('unimodaminoacids', 'UMAminoAcids'),
                      ('unimodmodifications', 'UMModifications'), ('unimodelements', 'UMElements'),
                      ('unimodspecificity', 'UMSpecificity'), ('config', 'Config'), ('bestprot', 'ParameterWTypes')]
        for t in tables:
            hdf.createTable(importGroup, t[0], t[1])

    def createETpeptidesTable(self):

        self.hdf.createTable(self.importGroup, 'etpeptides', 'ETpeptide')

    def writeParameters(self, paramDict, sectionName):
        """
        @brief writes parameter data using the provided section name
        @param paramDict <dictionary>: contaiing the parameter/value pairs
        @param sectionName <string>: section name for the parameters
        @return <dictionary>: giving success or error message
        """

        write = []
        for param in paramDict:
            write.append(dict(section=sectionName, parameter=param, value=paramDict[param]))

        return self.hdf.appendRows('/%s/parameters' % self.importGroup, write)

    def writeConfig(self, configList, table='config'):
        return self.hdf.appendRows('/%s/%s' % (self.importGroup, table), configList)

    def getMSconfig(self, cfg):
        readSections = ('xic', 'deisotoping')

        for rs in readSections:
            section = self.hdf.getDataEqual('/rawdata/config', 'set', rs)

            if rs in cfg.adjustedOptionalParams:
                params = {}
                for p in section:
                    if p['parameter'] not in cfg.adjustedOptionalParams[rs]:
                        params[p['parameter']] = p['value']
                if params:
                    cfg.parameters[rs].update(params)
            else:
                params = {}
                for p in section:
                    params[p['parameter']] = p['value']

                cfg.parameters[rs].update(params)

        cfg.convertParameters()

        return

    def getH5DFQuantMeth(self):
        try:
            isotopes = self.hdf.readTable('/rawdata/isotopes')
            if len(isotopes) == 0:
                return 0
            else:
                return isotopes[0]['method_id']
        except:
            return  0

    def writeMasses(self, massesDict):
        """
        @brief writes the masses data to hdf5 file
        @param massesDict <dictionary>: containing amino acid and other mass data
        @return <dictionary>: giving success or error message
        """

        write = []
        for name in massesDict:
            write.append(dict(name=name, mass=massesDict[name]))

        return self.hdf.appendRows('/%s/masses' % self.importGroup, write)

    def writeMods(self, modsDict):
        """
        @brief writes the Mascot modifications to the table
        @param modsDict <dictionary>: containing mascot amino acid definitions
        @return <dictionary>: giving success or error message
        """

        write = []
        for id in modsDict:
            data = modsDict[id]
            write.append(
                dict(id=id, name=data['name'], modtype=data['modtype'], da_delta=data['da_delta'], amino=data['amino'],
                     neutralloss=data['nl'], nlmaster=data['nlm'], relevant=data['relevant']))

        return self.hdf.appendRows('/%s/mods' % self.importGroup, write)

    def writeUnimodAminos(self, aaDict):
        """
        @brief writes the unimod amino acid dasta to the appropriate table
        @param aaDict <dictionary>: containing the amimo acid definitions
        @return <dictionary>: giving success or error message
        """

        write = []
        for aa in aaDict:
            data = aaDict[aa]
            compStr, compDict = self.calcComposition(data['elements'])
            write.append(dict(code=data['three_letter'], name=data['full_name'], title=data['title'],
                              monomass=data['mono_mass'], avgmass=data['avge_mass'], composition=compStr))

        return self.hdf.appendRows('/%s/unimodaminoacids' % self.importGroup, write)

    def writeUnimodElements(self, elements):
        """
        @brief writes the elemental data to the unimodelements table
        @param elements <dictionary>: containing elemental data
        @return <dictionary>: giving success or error message
        """

        write = []
        for e in elements:
            data = elements[e]
            write.append(dict(name=data['full_name'], title=data['title'],
                              monomass=data['mono_mass'], avgmass=data['avge_mass']))

        return self.hdf.appendRows('/%s/unimodelements' % self.importGroup, write)

    def getElementDictFromModsTable(self):
        modTable = self.hdf.tableTOC['/%s/unimodmodifications' % self.importGroup]
        elemDict = {}

        for c in modTable.colnames:
            if len(c) < 4:
                elemDict[c] = 0

        return elemDict

    def writeUnimodModifications(self, modifications):
        """
        @brief takes the unimod modification data and writes the basic data to the unimodmodifications table and
        the specificity data to the unimodspecificity table
        @param modifications <dictionary>: containing the unimod modifications
        @return <dictionary>: giving success or error message
        """

        writeMod = []
        writeSpec = []
        id = 0
        colsSpec = {}

        elementDict = self.getElementDictFromModsTable()
        colsMod = elementDict.copy()

        # generate the general data for modifications
        for mod in modifications:
            id += 1
            data = modifications[mod]
            name = unicodedata.normalize('NFKD', data['full_name']).encode('ascii', 'ignore')
            modDict = dict(modid=id, name=name, title=data['title'], delta_comp=data['delta']['composition'],
                           delta_mono_mass=data['delta']['mono_mass'], delta_avge_mass=data['delta']['avge_mass'])
            # modDict.update(elementDict)
            # modDict.update(compDict)
            # if 'Na' in modDict:
            #     del modDict['Na']
            # colsMod.update(modDict)
            writeMod.append(modDict)

            # generate the data for the modification specificity
            groups = data['specificity'].keys()
            groups.sort()
            for grp in groups:
                specify = data['specificity'][grp]
                specifyDict = dict(modid=id, group=grp, site=specify['site'], position=specify['position'],
                                   classification=specify['classification'], hidden=specify['hidden'])

                try:
                    for nl in specify['neutral_loss']:
                        specifyDict['nl_mono_mass'] = nl['mono_mass']
                        specifyDict['nl_avge_mass'] = nl['avge_mass']
                        specifyDict['nl_comp'] = nl['composition']
                except:
                    pass

                colsSpec.update(specifyDict)
                writeSpec.append(specifyDict)

        # write to hdf5
        self.hdf.appendRows('/%s/unimodmodifications' % self.importGroup, writeMod, colsMod.keys())
        self.hdf.appendRows('/%s/unimodspecificity' % self.importGroup, writeSpec)

        return

    def updateModsWithElements(self, uniMods):
        """
        @brief updates the mods table with the elemental composition from unimod
        @param uniMods <dictionary>: containing the unimod modification definitions
        @return <dictionary>: giving success or error message
        """

        modsTable = self.hdf.tableTOC['/%s/mods' % self.importGroup]

        for trow in modsTable:
            if trow['name'] in uniMods:
                uMod = uniMods[trow['name']]
                trow['composition'] = uMod['delta']['composition']
            else:
                raise ExHa.MascotModificationError('Modification "%s" not found in UniMod data' % trow['name'])
            trow.update()

    def calcComposition(self, elements):
        """
        @brief calculates a composition string and dictionary from the unumod elelments
        @param elements <dictionary>: containing the elemental composition
        @return <string>, <dictionary>:
        """
        keys = elements.keys()
        keys.sort()
        compStr = ''
        compDict = {}
        for k in keys:
            if elements[k] != 0:
                compStr += '%s ' % k
            else:
                compStr += '%s(%d) ' % (k, elements[k])

            if k == '13C':
                compDict['C13'] = elements[k]
            elif k == '15N':
                compDict['N15'] = elements[k]
            elif k == '2H':
                compDict['H2'] = elements[k]
            else:
                compDict[k] = elements[k]
        return compStr.strip(), compDict

    def writePeptides(self, peptides, pepTable):
        """
        @brief writes peptide and ETpeptide data to hdf5
        @param peptides <list>: containing peptide objects
        @return <dictionary>: giving success or error message
        """

        write = []

        for pep in peptides:
            # skip unwanted peptides
            if pep.retain == 0:
                continue

            pepDict = dict(query=pep.query, pepno=pep.pepno, sequence=pep.sequence, is_hook=pep.is_hook,
                           useinprot=pep.useinprot, modsVariable=pep.modsVariable, modsFixed=pep.modsFixed,
                           mass=pep.mass, da_delta=pep.da_delta, score=pep.score, misscleave=pep.misscleave,
                           numionsmatched=pep.numionsmatched, seriesfound=pep.seriesfound,
                           peaks1=pep.peaks1, peaks2=pep.peaks2, peaks3=pep.peaks3)

            if 'etMods' in pep.__dict__:
                pepDict['etModName'] = pep.etMods['desc']
                pepDict['etModDelta'] = pep.etMods['delta']

            write.append(pepDict)

        return self.hdf.appendRows('/%s/%s' % (self.importGroup, pepTable), write)

    def writeSeq2Acc(self, seq2acc):
        """
        @brief writes data tot he seq2acc table
        @param seq2acc <dictionary>: containing the protein data indexed by peptide sequence
        @return <dictionary>: giving success or error message
        """

        write = []

        for seq in seq2acc:
            # get peptide data first
            data = seq2acc[seq]
            pepDict = dict(sequence=seq, hook=data['hook'], bestczrank=data['bestczrank'],
                           pepscore=data['pepscore'], hookscore=data['hookscore'], numpeps=data['numpep'])

            # now loop through the protein accessions and add a line for each accession
            for acc in data['prots']:
                hittype = 'FWD'
                if acc['accession'].find('###REV###') >-1 or acc['accession'].find('###RND###') >-1 \
                        or acc['accession'].startswith('DD'):
                    hittype = 'REV'
                protDict = dict(accession=acc['accession'], start=acc['start'], end=acc['end'], hittype=hittype)
                protDict.update(pepDict)
                write.append(protDict)

        return self.hdf.appendRows('/%s/seq2acc' % self.importGroup, write)

    def getPeptidesFromProtein(self, proteinAcc):
        """
        @brief reads the seq2acc table to find all peptides linked to the protein accession
        """

        # fetch the required data
        pepArray = self.hdf.getDataEqual('/%s/seq2acc' % self.importGroup, 'accession', proteinAcc)
        headers = self.hdf.readTable('/rawdata/msmsheader')
        spec2head = {}
        for idx, h in enumerate(headers):
            spec2head[h['spec_id']] = idx
        tmp = self.hdf.readTable('/%s/queries' % self.importGroup)
        query2spec = {}
        for t in tmp:
            query2spec[t['query']] = t['spec_id']
        tmp = 0

        peptides = {}
        fwhmList = []
        rtList = []

        for pep in pepArray:
            clause = '(sequence == "%s") & (useinprot == 1)' % pep['sequence']
            foundPeps = self.hdf.getDataGeneral('/%s/peptides' % self.importGroup, clause)

            peptideList = []
            for fp in foundPeps:
                spec = query2spec[fp['query']]
                head = headers[spec2head[spec]]
                peptideList.append(dict(query=fp['query'], rank=fp['pepno'], spec_id=spec, ishook=fp['is_hook'],
                                        score=fp['score'], fwhm=head['fwhm'], rt=head['rtapex']))

            peptides[pep['sequence']] = dict(hssHook=pep['hook'], hookScore=pep['hookscore'], pepScore=pep['pepscore'],
                                             bestRank=pep['bestczrank'], numPeps=pep['numpeps'], peptides=peptideList)

        return peptides

    def getPeptidesFromProteinNoMSMSheader(self, proteinAcc):
        """
        @brief reads the seq2acc table to find all peptides linked to the protein accession
        """

        # fetch the required data
        pepArray = self.hdf.getDataEqual('/%s/seq2acc' % self.importGroup, 'accession', proteinAcc)
        tmp = 0

        peptides = {}
        fwhmList = []
        rtList = []

        for pep in pepArray:
            clause = '(sequence == "%s") & (useinprot == 1)' % pep['sequence']
            foundPeps = self.hdf.getDataGeneral('/%s/peptides' % self.importGroup, clause)

            peptideList = []
            for fp in foundPeps:
                peptideList.append(dict(query=fp['query'], rank=fp['pepno'], ishook=fp['ishook'], score=fp['score']))

            peptides[pep['sequence']] = dict(hssHook=pep['hook'], hookScore=pep['hookscore'], pepScore=pep['pepscore'],
                                             bestRank=pep['bestczrank'], numPeps=pep['numpeps'], peptides=peptideList)

        return peptides

    def writeStatistics(self, statistics):
        """
        @brief writes the summary statistics to hdf5
        @param statistics <dictionary>: containing statistic value pairs
        @return <dictionary>: giving success or error message
        """

        write = []
        for stat in statistics:
            write.append(dict(statistic=stat, value=float(statistics[stat])))

        return self.hdf.appendRows('/%s/statistics' % self.importGroup, write)

    def getNameFromProteinAcc(self, proteinAcc):

        prot = self.hdf.getDataEqual('/%s/proteins' % self.importGroup, 'accession', proteinAcc)
        if len(prot) == 0:
            return ''
        else:
            return prot[0]['name']

    def writeBestProtein(self, bestProt):
        """
        @brief writes the protein data for the best protein
        @param bestProt <dictionary>: containing the protein data
        """

        dataList = []
        for key in bestProt:
            if isinstance(bestProt[key], float):
                dataType = 'float'
            elif isinstance(bestProt[key], int):
                dataType = 'int'
            else:
                dataType = 'str'
            dataList.append(dict(parameter=key, type=dataType, value=str(bestProt[key])))

        return self.hdf.appendRows('/%s/bestprot' % self.importGroup, dataList)

    def readBestProtein2Dict(self):
        bestProtDict = {}
        # get the summary data
        data = self.hdf.readTable('/%s/bestprot' % self.importGroup)
        for d in data:
            if d['type'] == 'int':
                bestProtDict[d['parameter']] = int(d['value'])
            elif d['type'] == 'float':
                bestProtDict[d['parameter']] = float(d['value'])
            else:
                bestProtDict[d['parameter']] = d['value']

        return bestProtDict

    def writeProteins(self, protList):
        """
        @brief writes the list of proteins directly to the hdf5 file
        @param protList <list>: containing protien dictionaries
        @return <dictionary>: giving success or error message
        """

        return self.hdf.appendRows('/%s/proteins' % self.importGroup, protList)

    def writeIndex(self, indexList):

        write = []
        for item in indexList:
            write.append(dict(section=item[0], linenumber=item[1]))

        return self.hdf.appendRows('/%s/index' % self.importGroup, write)

    def writeQueries(self, queries):

        write = []
        for id in queries:
            data = queries[id]
            queryDict = dict(query=id, msms_id=data['msmsid'], spec_id=int(data['msmsid'][1:]), rt=data['rt'],
                             prec_neutmass=data['prec_neutmass'], prec_mz=data['prec_mz'],
                             prec_charge=data['prec_charge'], homology=data['homology'], matches=data['matches'],
                             numpeps=data['numpeps'], delta_seq=0.0, delta_mod=0.0)
            if data['numpeps'] > 0:
                queryDict['delta_seq'] = data['delta_seq']
                queryDict['delta_mod'] = data['delta_mod']
            write.append(queryDict)

        return self.hdf.appendRows('/%s/queries' % self.importGroup, write)

    def writeCalMasses(self, data):
        write = []
        queryList = data.keys()
        queryList.sort()

        names = self.hdf.readTable('/calibration/calmasses').dtype.names

        for query in queryList:
            if 'pepno' in data[query]:
                calDict = {}
                for name in names:
                    calDict[name] = data[query][name]
                write.append(calDict.copy())
        return self.hdf.appendRows('/calibration/calmasses', write)

    def createTableIndexes(self, isPreCal):
        # update imports status
        importsTable = self.hdf.tableTOC['/imports']
        importsTable.cols.status[self.impRow] = 'completed'
        ig = self.importGroup

        if isPreCal:
            indexes = [('/%s/calmasses' % ig, ['query', 'sequence', 'spec_id'])]
        else:
            indexes = [('/%s/queries' % ig, ['query']),
                       ('/%s/peptides' % ig, ['query', 'sequence', 'is_hook', 'pepno', 'score']),
                       ('/%s/proteins' % ig, ['accession']), ('/%s/seq2acc' % ig, ['accession', 'sequence'])]

        # create indexes
        for index in indexes:
            self.hdf.indexTable(index[0], index[1])

        return

    def readModDicts(self):
        """
        @brief returns all the modifications
        """

        data = self.hdf.readTable('/' + self.importGroup + '/mods')
        fixed = {}
        variable = {}
        for d in data:
            if 'composition' in data.dtype.names:
                # has composition field so calculate compostion from this
                elementsDict = self.compostition2dict(d['composition'])
            else:
                elementsDict = {'C': d['C'], 'H': d['H'], 'N': d['N'], 'O': d['O'], 'P': d['P'],
                                '13C': d['C13'], '2H': d['H2'], '15N': d['N15']}

            mod = dict(name=d['name'], delta=d['da_delta'], amino=d['amino'], relevant=d['relevant'],
                       nlmaster=d['nlmaster'], neutralloss=d['neutralloss'], elem=elementsDict)

            if d['modtype'] == 'fixed':
                fixed[d['id']] = mod.copy()
            elif d['modtype'] == 'variable':
                variable[d['id']] = mod.copy()

        return variable, fixed

    def readIsotopes(self):

        data = self.hdf.readTable('/rawdata/isotopes')
        isotopeList = []
        for d in data:
            isotope = {}
            for c in d.dtype.names:
                isotope[c] = d[c]

            isotopeList.append(isotope.copy())

        return isotopeList

    def compostition2dict(self, composition):

        bits = composition.split(' ')
        elements = {}
        for element in bits:
            if '(' in element:
                elem, num = element[:-1].split('(')
                elements[elem] = int(num)
            else:
                elements[element] = 1

        return elements

    def getMasses(self):
        hdf = self.hdf
        tablePath = '/%s/%s' % (self.importGroup, 'masses')
        masses = hdf.readTable(tablePath)
        return masses
