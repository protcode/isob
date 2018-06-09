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

This file is the primary handler for hdf5 files.  Dealing with the creation
of groups, tables & indexes and the entry, retrieval and deletion of data.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob

"""
import numpy as np
from numpy.lib import recfunctions as rfn
import tables as tables
from tables.nodes import filenode
from pathlib import Path
import hdf5Tables
import ExceptionHandler as ExHa


class hdf5Base:
    """
    @brief class to read and write pyMSsafe data into HDF5 format
    """

    def __init__(self, hdf5file, overwrite=False, forceCreation=False):
        self.filePath = Path(hdf5file)
        self.filePathStr = str(hdf5file)
        self.filt = tables.Filters(complevel=5)
        self.isopen = 0

        # if file exists read TOC data and if not existing create new file
        exists = self.filePath.exists()
        if exists and not overwrite:
            self.openHDF('r')
            self.close()
        elif (exists and overwrite) or (not exists and forceCreation):
            self.tableTOC = {}
            self.groupTOC = {}
            self.arrayTOC = {}
            self.hdf5 = tables.openFile(self.filePathStr, mode='w', filters=self.filt)
            self.hdf5.close()
        elif not exists:
            raise ExHa.HDFwritingError('HDF5 does not exist:', str(self.filePath))

    def loadTableTOC(self):
        allTables = {}
        for table in self.hdf5.walkNodes(classname='Table'):
            # noinspection PyProtectedMember
            allTables[table._v_pathname] = table
        self.tableTOC = allTables
        return allTables

    def loadEArrayTOC(self):
        allEArrays = {}
        for eArray in self.hdf5.walkNodes(classname='EArray'):
            # noinspection PyProtectedMember
            allEArrays[eArray._v_pathname] = eArray
        self.arrayTOC = allEArrays
        return allEArrays

    def loadGroupTOC(self):
        allGroups = {}
        for group in self.hdf5.walkGroups():
            # noinspection PyProtectedMember
            allGroups[group._v_name] = group

        self.groupTOC = allGroups
        return allGroups

    def openHDF(self, mode):
        """
        @brief opens an existing hdf5 file
        """

        self.hdf5 = tables.openFile(self.filePathStr, mode=mode, filters=self.filt)

        self.groupTOC = self.loadGroupTOC()
        self.arrayTOC = self.loadEArrayTOC()
        self.tableTOC = self.loadTableTOC()

        self.isopen = mode

    def readOpen(self):
        """
        @brief opens an existing hdf5 file
        """
        self.openHDF('r')

    def appendOpen(self):
        """
        @brief opens an existing hdf5 file
        """
        self.openHDF('a')

    def close(self):
        self.hdf5.close()
        self.isopen = 0

    @staticmethod
    def joinGroup2Table(group, table):
        if group == '/':
            tablePath = '/%s' % table
        else:
            tablePath = '/%s/%s' % (group, table)
        return tablePath

    def createGroup(self, group):
        """
        @breif creates a new group.  if the group already exists then nothing happens
        @param group <string>: the name of the group the table belongs to
        """
        if group in self.groupTOC:
            # group exists
            raise ExHa.HDF5consistancyError('Group %s already exists' % group)
        else:
            newGroup = self.hdf5.createGroup(self.hdf5.root, group)
            # noinspection PyProtectedMember
            self.groupTOC[newGroup._v_name] = newGroup

    def removeGroup(self, group):
        """
        @breif removes a group from the hdf5 file
        @param group <string>: the name of the group the table belongs to
        """

        if group in self.groupTOC:
            # group exists
            groupObj = self.groupTOC.pop(group)
            # noinspection PyProtectedMember
            groupObj._f_remove(recursive=True)

            # remove entries for deleted tables
            remTables = [x for x in self.tableTOC if group in x]
            for rt in remTables:
                del self.tableTOC[rt]

    def createTable(self, group, name, definition, delete=False):
        """
        @breif creates a new table.  if the table already exists then optionally the old version can be deleted
        @param group <string>: the name of the group the table belongs to
        @param name <string>: the name of the new table
        @param definition <string>: the definition class of the table
        @param delete <boolean>: flag to control the deletion of old tables
        """
        tablePath = self.joinGroup2Table(group, name)

        # First find if the table already exists
        if tablePath in self.tableTOC:
            # table exists, do we delete it
            if delete:
                # we can delete it
                # noinspection PyProtectedMember
                self.tableTOC[tablePath]._f_remove()
            else:
                # return an error
                raise ExHa.HDF5consistancyError('Table %s exists, not deleted' % tablePath)

        if hasattr(hdf5Tables, definition):
            tClass = getattr(hdf5Tables, definition)
        else:
            raise ExHa.HDF5consistancyError('Table definition "%s" does not exist' % definition)

        if group in self.groupTOC:
            # group is present so create the table
            newTable = self.hdf5.createTable(self.groupTOC[group], name, tClass)
            # noinspection PyProtectedMember
            self.tableTOC[newTable._v_pathname] = newTable
            if hasattr(hdf5Tables, definition + '_attr'):
                newAttr = getattr(hdf5Tables, definition + '_attr')
            else:
                newAttr = {}
            newTable.attrs.columnDesc = newAttr

        else:
            raise ExHa.HDF5consistancyError('Group %s does not exist, table not created' % group)

    def getColumnDescAttr(self, tablePath):
        """
        @brief gets the attribute columnDesc from the table
        @param tablePath <string>: the path name of the new table
        """
        if tablePath in self.tableTOC:
            tableObj = self.tableTOC[tablePath]
            if hasattr(tableObj.attrs, 'columnDesc'):
                return tableObj.attrs.columnDesc
            else:
                return []

    def getColumnDescDict(self, tablePath):
        """
        @brief gets the attribute columnDesc from the table
        @param tablePath <string>: the path name of the new table
        """
        if tablePath in self.tableTOC:
            tableObj = self.tableTOC[tablePath]
            if hasattr(tableObj.attrs, 'columnDesc'):
                cdDict = {}

                attList = tableObj.attrs.columnDesc
                for data in attList:
                    col, desc = data[4:].split(': ', 1)
                    cdDict[col] = desc

                for cn in tableObj.dtype.descr:
                    if cn[0] not in cdDict:
                        # missing column description
                        cdDict[cn[0]] = '-'
                return cdDict
            else:
                return {}

    def getColumnDescListWithType(self, tablePath):
        """
        @brief gets the attribute columnDesc from the table
        @param tablePath <string>: the path name of the new table
        """
        oList = []
        if tablePath in self.tableTOC:
            tableObj = self.tableTOC[tablePath]
        attDict = self.getColumnDescDict(tablePath)

        for d in tableObj.dtype.descr:
            if d[1].startswith('<f'):
                # float type
                bytes = int(d[1][2:])
                dt = 'float %ibit' % (bytes * 8)
            elif d[1].startswith('<i'):
                # integer type
                bytes = int(d[1][2:])
                dt = 'integer %ibit' % (bytes * 8)
            elif d[1].startswith('|S'):
                dt = 'string %s characters' % d[1][2:]
            else:
                # print 'missing dtype', d[1]
                dt = 'unknown'
            oList.append((d[0], dt, attDict[d[0]]))
        return oList

    def removeTable(self, tablePath):
        """
        @brief deletes an existing table table.  if the table already exists then the old version will be deleted
        @param tablePath <string>: the path name of the new table
        """

        # First find if the table already exists
        if tablePath in self.tableTOC:
            # table exists, we delete it
            # noinspection PyProtectedMember
            self.tableTOC[tablePath]._f_remove()
            del self.tableTOC[tablePath]

    def renameTable(self, group, oldName, newName):
        """
        @brief renames an existing table table.  if the table already exists then the old version will be deleted
        @param tablePath <string>: the path name of the new table
        """
        oldTablePath = '/%s/%s' % (group, oldName)
        newTablePath = '/%s/%s' % (group, newName)
        # First find if the table already exists
        if oldTablePath in self.tableTOC:
            # table exists, we rename it
            table = self.tableTOC[oldTablePath]
            table._f_rename(newName)

            # update the tableTOC entries
            del self.tableTOC[oldTablePath]
            self.tableTOC[newTablePath] = table

    def indexTable(self, tablePath, columnsList):

        # test that the table is present
        if tablePath not in self.tableTOC:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)
        else:
            # create the table object for data appending
            table = self.tableTOC[tablePath]

        # flush the table data then create the indexes
        table.flush()
        for col in columnsList:
            if col in table.colnames:
                colInstance = table.colinstances[col]
                colInstance.createIndex()
                table.flush()
            else:
                raise ExHa.HDF5consistancyError('Table %s does not have column %s' % (tablePath, col))

        return dict(code=0)

    def appendRows(self, tablePath, dataList, keyList=None):
        """
        @brief appends data for a single row to the table
        @param tablePath <string>: internal path of the table to be read.
        @param dataList <list>: containing dictionaries holding required data (NB also can be a single ndArray)
        @param keyList <list>: optionally containing a list of keys used in the data: use for dataLists with variable
            keys for each data dictionary
        """
        if not keyList:
            keyList = []

        # test if any data to be written
        if len(dataList) == 0:
            raise ExHa.HDFmissingDataError('No data to be written to %s' % tablePath)

        # test that the table is present
        if tablePath not in self.tableTOC:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)
        else:
            # create the row object for data appending
            table = self.tableTOC[tablePath]
            rowObj = table.row

        # test the fields in the dataList for compatibility with the table fields
        useDict = 1
        if keyList:
            keys = keyList
            # useDict = 0
        else:
            test = dataList[0]
            if isinstance(test, dict):
                keys = test.keys()
            else:
                keys = test.dtype.names
                useDict = 0

        for key in keys:
            if key not in table.colnames:
                raise ExHa.HDF5consistancyError('Table %s does not have column %s' % (tablePath, key))

        # append data row by row
        if useDict:
            # dictionary can itarate over the keys
            for data in dataList:
                keyList = [x for x in data if x in keys]
                for key in keyList:
                    rowObj[key] = data[key]
                rowObj.append()
        else:
            for data in dataList:
                # ndarray itarates over the data so must supply the keys
                for key in keys:
                    rowObj[key] = data[key]
                rowObj.append()

        # flush the buffer
        table.flush()
        return

    def completely_sorted_index_table(self, table_path, col):
        # test that the table is present
        if table_path not in self.tableTOC:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % table_path)
        else:
            # create the table object for data appending
            table = self.tableTOC[table_path]

        # flush the table data then create the indexes
        table.flush()
        if col in table.colnames:
            col_instance = table.colinstances[col]
            # print 'creating sorted index on %s' % col

            col_instance.createCSIndex()
            table.flush()
        else:
            raise ExHa.HDF5consistancyError('Table %s does not have column %s' % (table_path, col))
        return dict(code=0)

    def updateDataInRow(self, tablePath, whereClause, changeData):
        """
        @brief changes specific data in a table.  All rows matching rowFilter will be changed
        @param tablePath <string>: internal path of the table to be read.
        @param whereClause <string>: filter conditions for the table.
        @param changeData <dictionary>: containing new data values indexed by the column name
        """

        # test that the table is present
        if tablePath not in self.tableTOC:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)
        else:
            # create the row object for data appending
            table = self.tableTOC[tablePath]

            for row in table.where(whereClause):
                for col in changeData:
                    row[col] = changeData[col]
                row.update()

            table.flush()

    def removeRows(self, tablePath, conditions):
        """
        @brief appends data for a single row to the table
        @parmam tablePath <string>: internal path of the table to be read.
        @parmam conditions <dictionary>: containing filter values indexed by the column name
        """

        # test that the table is present
        if tablePath not in self.tableTOC:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)
        else:
            # create the row object for data appending
            table = self.tableTOC[tablePath]

            rows2remove = []

            for idx in xrange(table.nrows):
                remove = True
                for key in conditions:
                    if not table[idx][key] == conditions[key]:
                        remove = False
                        break
                if remove:
                    rows2remove.append(idx)

            rows2remove.sort(reverse=True)

            for idx in rows2remove:
                table.removeRows(idx)

    def hasTable(self, tablePath):
        """
        @brief returns True if hdf5 file contains a table at the specified path, otherwise returns False
        @parmam tablePath <string>: internal path of the table to be checked.
        """
        return tablePath in self.tableTOC

    def getNumberOfTableRows(self, tablePath):
        """
        @brief returns True if hdf5 file contains a table at the specified path, otherwise returns False
        @parmam tablePath <string>: internal path of the table to be checked.
        """
        if tablePath in self.tableTOC:
            return self.tableTOC[tablePath].nrows
        else:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)

    def getTable(self, tablePath):
        """
        @brief returns all the data from hdf5 table specified by tablePath
        @parmam tablePath <string>: internal path of the table to be read.
        @return data <array>: containing the table contents
        """
        if tablePath in self.tableTOC:
            table = self.tableTOC[tablePath]
            return table
        else:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)

    def readTable(self, tablePath):
        """
        @brief returns all the data from hdf5 table specified by tablePath
        @parmam tablePath <string>: internal path of the table to be read.
        @return data <array>: containing the table contents
        """
        if tablePath in self.tableTOC:
            table = self.tableTOC[tablePath]
            return table.read()
        else:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)

    def getDataInKernelSearch(self, tablePath, whereClause, column):
        """
        @brief performs an in-kernel search of tablePath for data matching the whereClause, returns
        non redundant set of single value data given by column. It is beneficial since the search is
        in-kernel and not all the data is read into memory
        @return set: set of contents of column
        @param tablePath <string>: internal path of the table to be read.
        @param whereClause <string>: filter conditions for the table.
        @param column <string>: single column value from table selected by whereClause.
        @return data <array>: containing the table contents
        """

        if tablePath in self.tableTOC:
            table = self.tableTOC[tablePath]
            return set([x[column] for x in table.where(whereClause)])
        else:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)

    def getDataGeneral(self, tablePath, whereClause):
        """
        @brief searches tablePath for data matching the whereClause, returns all matching rows
        @return data <ndarray>: containing the msmsheader data
        @parmam tablePath <string>: internal path of the table to be read.
        @parmam whereClause <string>: filter conditions for the table.
        @return data <array>: containing the table contents
        """

        if tablePath in self.tableTOC:
            table = self.tableTOC[tablePath]
            return table.readWhere(whereClause)
        else:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)

    def getDataEqual(self, tablePath, column, value):
        """
        @brief searches column in tablePath for data matching value, returns all matching rows
        @return data <ndarray>: containing the msmsheader data
        @parmam tablePath <string>: internal path of the table to be read.
        @parmam column <string>: filter conditions for the table.
        @parmam value <any>: filter conditions for the table, should be string or integer or float.
        @return rtn <array>: containing the table contents
        """
        if isinstance(value, basestring):
            clause = '(%s == "%s")' % (column, value)
        else:
            clause = '(%s == %s)' % (column, str(value))

        rtn = self.getDataGeneral(tablePath, clause)

        return rtn

    def getDataRange(self, tablePath, column, low, high):
        """
        @brief searches column in tablePath for data matching value, returns all matching rows
        @return data <ndarray>: containing the msmsheader data
        @parmam tablePath <string>: internal path of the table to be read.
        @parmam column <string>: filter conditions for the table.
        @parmam value <string>: filter conditions for the table.
        @parmam value <string>: filter conditions for the table.
        @return rtn <array>: containing the table contents
        """
        clause = '(%(col)s >= %(low)s) & (%(col)s <= %(high)s)' % {'col': column, 'low': str(low), 'high': str(high)}

        rtn = self.getDataGeneral(tablePath, clause)

        return rtn

    def replaceTable(self, group, table, newData, definition):
        """
        @brief deletes an existing table and replaces it with a new table containing newData
        """
        tablePath = self.joinGroup2Table(group, table)

        if tablePath in self.tableTOC:
            tableObj = self.tableTOC[tablePath]
            indexesDict = tableObj.colindexed
            indexes = [x for x in indexesDict if indexesDict[x]]

            self.createTable(group, table, definition, True)
            self.appendRows(tablePath, newData)

            self.indexTable(tablePath, indexes)
        else:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)

    def createNodeWithFileData(self, group, nodeName, fileHandle):
        fNode = filenode.newNode(self.hdf5, where=group, name=nodeName)
        fNode.write(fileHandle.read())
        fNode.close()
        self.hdf5.flush()

        self.loadEArrayTOC()

    def removeNode(self, nodePath):
        if nodePath in self.arrayTOC:
            # node exists, we delete it
            # noinspection PyProtectedMember
            self.arrayTOC[nodePath]._f_remove()
            del self.arrayTOC[nodePath]

    def createFileNodeForImage(self, group, nodeName, fileHandle):
        # find any fileNodes in the current group
        nodePath = self.joinGroup2Table(group, nodeName)
        whereGroup = '/%s' % group
        if nodePath in self.arrayTOC:
            self.hdf5.removeNode(self.arrayTOC[nodePath])

        self.createNodeWithFileData(whereGroup, nodeName, fileHandle)
        self.loadEArrayTOC()

        self.hdf5.flush()

    def getNodes(self, nodeNamesToGet, className, group):
        """get all nodes of a certain type matching names in nodenamestoget list"""

        importNode = self.groupTOC[group]

        # noinspection PyProtectedMember
        nodeList = importNode._f_listNodes(className)
        fetchedNodeList = [x for x in nodeList if x.name in nodeNamesToGet]

        return fetchedNodeList

    def table2file(self, tablePath, fileDir=False):
        """
        @brief takes the contents of a table and outputs the data to a file of the same name
        @parmam tablePath <string>: internal path of the table to be read.
        """

        if tablePath in self.tableTOC:
            table = self.tableTOC[tablePath]
            data = table.read()
            if fileDir:
                filePath = Path(fileDir, Path(tablePath).name + '.txt')
            else:
                filePath = self.filePath.parent.joinpath(Path(tablePath).name + '.txt')
            titles = data.dtype.names
            fout = open(str(filePath), 'w')

            fout.write('\t'.join(titles) + '\n')

            for d in data:
                lineList = []
                for name in titles:
                    lineList.append(str(d[name]))
                fout.write('\t'.join(lineList) + '\n')

            fout.close()
            return filePath
        else:
            raise ExHa.HDF5consistancyError('Table %s does not exist' % tablePath)

    def checkTableAgainstDefinition(self, inTable, definition, adapt=False):
        """
        @brief Checks wether the columns in the numpy array 'inTable' suit those in the table definition from
                CommonUtils/hdf5Tables. If the 'adapt' flag is set to True, missing columns will be added to the table
                and extra columns will be removed and the adjusted table is returned. Otherwise the function returns a
                logical value indicating whether the given numpy array suits the definition (True) or not (False).
        @param inTable: <ndarray> numpy array to be checked and adapted if necessary
        @param definition: <String> name of the table definition CommonUtils/hdf5Tables to check against
        @param adapt: <logical> flag indicating whether the given numpy array shall be adapted to the definition and
                the modified numpy array shal be returned (default=False)
        """

        suitsDefinition = True

        if hasattr(hdf5Tables, definition):
            tmpTabObj = hdf5Tables.__dict__[definition]
        else:
            raise ExHa.configError("Can't find table %s in CommonUtils/hdf5Tables", definition)

        definedCols = tmpTabObj.columns.keys()
        colsInTable = list(inTable.dtype.names)

        extraColsInTable = list(set(colsInTable) - set(definedCols))
        missingColsInTable = list(set(definedCols) - set(colsInTable))
        intersect = list(set(definedCols) & set(colsInTable))

        if len(extraColsInTable) > 0:
            suitsDefinition = False
            print "Given table '%s' has the following extra columns: %s" % (definition, ', '.join(extraColsInTable))
            if adapt:
                print "  -> adapting table '%s'; removing extra columns. " % definition
                inTable = inTable[intersect]

        if len(missingColsInTable) > 0:
            suitsDefinition = False
            print "Column(s) %s mssing in the given table '%s'" % (', '.join(missingColsInTable), definition)
            if adapt:
                for mcit in missingColsInTable:
                    print "  -> adding empty column '%s'" % mcit
                    inTable = rfn.append_fields(inTable, names=mcit,
                                                data=np.asarray([None] * inTable.size,
                                                                dtype=tmpTabObj.columns[mcit].dtype),
                                                usemask=False)
        if adapt:
            return inTable
        else:
            return suitsDefinition
