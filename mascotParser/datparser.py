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

This file mainly handles the reading of Mascot .dat files grouping the data
by the sections defined by Mascot.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""

import sys
import re

# CommonUtils imports
sys.path.insert(0, '..')
import CommonUtils.ExceptionHandler as ExHa

# project imports

rx_sectionh = re.compile('name="(?P<name>\w+)"')
rx_paramline = re.compile('(?P<key>^[A-Za-z0-9_]+)=(?P<value>.*)')
rx_peptide = re.compile('(?P<key>^q[0-9]+_p[0-9]+(_[a-z_]*)?)=(?P<value>.*)')
rx_query = re.compile('query(?P<ID>[0-9]+)')
rx_summary = re.compile('(?P<key>(qmass|qexp|qmatch|qplughole)[0-9]*)=(?P<value>.*)')
rx_enzymeline = re.compile('(?P<key>^[A-Za-z0-9_]+):(?P<value>.*)')
rx_protein = re.compile('"(?P<key>.+)"=(?P<value>.*)')


class DatParser:
    def __init__(self, datfile, cfg, logs):
        """
        __init__()

        initialises class and creates links to the logger and dat file oject
        """
        self.datfileobj = datfile
        self.cfg = cfg
        self.logs = logs
        self.sectionsNeeded = {}

    def startParsing(self):
        """
        startParsing(string) -> integer

        Parses the dat file given by datfilepath.
        """
        datPath = str(self.datfileobj.filedic['datpath'])
        self.logs.datlog.info('Starting parsing file %s' % datPath)
        fin = open(datPath, 'r')

        try:
            self.cargo = fin.next()
        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Empty dat file.')
            raise

        statedic = {'start':          self.doSkip(fin, self.datfileobj),
                    'parameters':     self.doParameters(fin, self.datfileobj),
                    'masses':         self.doMasses(fin, self.datfileobj),
                    'unimod':         self.doUnimod(fin, self.datfileobj),
                    'enzyme':         self.doEnzyme(fin, self.datfileobj),
                    'header':         self.doHeader(fin, self.datfileobj),
                    'summary':        self.doSummary(fin, self.datfileobj),
                    'et_summary':     self.doETsummary(fin, self.datfileobj),
                    'peptides':       self.doPeptides(fin, self.datfileobj),
                    'postpeptides':   self.doPostPeptides(fin, self.datfileobj),
                    'et_peptides':    self.doETpeptides(fin, self.datfileobj),
                    'proteins':       self.doProteins(fin, self.datfileobj),
                    'query':          self.doQuery(fin, self.datfileobj),
                    'postquery':      self.doSkip(fin, self.datfileobj),
                    'taxonomy':       self.doSkip(fin, self.datfileobj),
                    'index':          self.doIndex(fin, self.datfileobj),
                    'quantitation':   self.doSkip(fin, self.datfileobj),
                    'boundary':       self.doBoundary(fin, self.datfileobj)}

        for state in statedic:
            if statedic[state].gi_code.co_name in ('doSkip', 'doBoundary'):
                pass
            elif state.startswith('et_'):
                pass
            else:
                self.sectionsNeeded[state] = 0

        coroutine = 'start'
        try:
            while 1:
                (coroutine, self.cargo) = statedic[coroutine].next()
        except StopIteration:
            self.doPostParsing(self.datfileobj)
            self.logs.datlog.info('Parsing complete')
        except:
            raise
        fin.close()

    def isBoundary(self, text):
        """
        @brief checks for the boundary condition (always preceeded by comment marks)
        """
        return text[:2] == '--'

    def doSkip(self, fileiterator, datfile):
        """
        @brief reads through the file until the next boundry
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        isbound = self.isBoundary

        while 1:
            while not isbound(self.cargo):
                self.cargo = fileiterator.next()
            self.cargo = fileiterator.next()
            yield ('boundary', self.cargo)

    def doParameters(self, fileiterator, datfile):
        """
        @brief reads the parameters section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing parameters')
            while 1:
                allparams = {}
                while not isbound(self.cargo):
                    match = rx_paramline.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')
                    self.cargo = fileiterator.next()
                datfile.hdfMascot.writeParameters(allparams, 'parameters')
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doParameters')
            raise

    def doMasses(self, fileiterator, datfile):
        """
        @brief reads the masses section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """

        try:
            isbound = self.isBoundary
            self.logs.datlog.info('Loading and Parsing masses')
            while 1:
                allparams = {}
                while not isbound(self.cargo):
                    match = rx_paramline.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')
                    self.cargo = fileiterator.next()
                datfile.addMasses(allparams)
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doMasses')
            raise

    def doUnimod(self, fileiterator, datfile):
        """
        @brief reads the masses section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """

        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing Unimod data')
            while 1:
                xml = []
                while not isbound(self.cargo):
                    xml.append(self.cargo)
                    self.cargo = fileiterator.next()
                datfile.addUnimod(xml)
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doUnimod')
            raise

    def doEnzyme(self, fileiterator, datfile):
        """
        @brief reads the enzyme section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing Enzyme')
            while 1:
                allparams = {}
                while not isbound(self.cargo):
                    match = rx_enzymeline.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')
                    elif len(self.cargo) > 2:
                        allparams['Side'] = self.cargo[:-1]
                    self.cargo = fileiterator.next()
                datfile.hdfMascot.writeParameters(allparams, 'enzyme')
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doEnzyme')
            raise

    def doHeader(self, fileiterator, datfile):
        """
        @brief reads the header section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            isbound = self.isBoundary
            self.logs.datlog.info('Loading and Parsing header')
            while 1:
                allparams = {}
                while not isbound(self.cargo):
                    match = rx_paramline.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')
                    self.cargo = fileiterator.next()
                datfile.hdfMascot.writeParameters(allparams, 'header')
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doHeader')
            raise

    def doSummary(self, fileiterator, datfile):
        """
        @brief reads the summary section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing summary')
            while 1:
                allparams = {}
                while not isbound(self.cargo):
                    match = rx_summary.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')

                    self.cargo = fileiterator.next()
                datfile.addSummary(allparams)
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doSummary')
            raise

    def doETsummary(self, fileiterator, datfile):
        """
        @brief reads the summary section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing et_summary')
            while 1:
                allparams = {}
                while not isbound(self.cargo):
                    match = rx_summary.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')

                    self.cargo = fileiterator.next()
                datfile.addETsummary(allparams)
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doETsummary')
            raise

    def doPeptides(self, fileiterator, datfile):
        """
        @brief reads the peptide section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        key = 'None'
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing peptides')
            collection = {}
            lastquery = ''
            while 1:
                while not isbound(self.cargo):
                    match = rx_peptide.search(self.cargo)
                    if match:
                        value = match.group('value')
                        key = match.group('key')
                        if value == '-1':
                            # no matching peptide
                            datfile.stats['numspectra_nopeps'] += 1
                        else:
                            query = key.split('_')[0]
                            if query != lastquery:
                                lastquery = query
                                if collection:
                                    # add the query number to addpeptides
                                    datfile.addPeptides(collection)
                                    collection = {}
                            collection[key] = value
                    self.cargo = fileiterator.next()
                    # add last peptide
                if collection:
                    # add the query number to addpeptides
                    datfile.addPeptides(collection)
                    pass
                yield ('postpeptides', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'doPeptides, last peptide = %s' % key)
            raise

    def doETpeptides(self, fileiterator, datfile):
        """
        @brief reads the peptide section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        key = 'None'
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing et_peptides')
            collection = {}
            lastquery = ''

            # create the etpeptides table in the HDF5 file
            datfile.createETpeptidesTable()
            while 1:
                while not isbound(self.cargo):
                    match = rx_peptide.search(self.cargo)
                    if match:
                        value = match.group('value')
                        key = match.group('key')
                        if value.startswith('-1'):
                            # no matching peptide
                            datfile.stats['numspectra_nopeps'] += 1
                        else:
                            query = key.split('_')[0]
                            if query != lastquery:
                                lastquery = query
                                if collection:
                                    # add the query number to addpeptides
                                    datfile.addETpeptides(collection)
                                    collection = {}
                            collection[key] = value
                    self.cargo = fileiterator.next()
                    # add last peptide
                if collection:
                    # add the query number to addpeptides
                    datfile.addETpeptides(collection)
                yield ('postpeptides', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'doETpeptides, last peptide = %s' % key)
            raise

    def doPostPeptides(self, fileiterator, datfile):
        """
        @brief processes the global peptide data
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        while 1:
            self.logs.datlog.info('%d peptides processed' % datfile.peptidecounter)
            self.logs.datlog.info('%d failed peptides' % datfile.failedsequences)
            datfile.doPostPeptides()
            yield ('boundary', self.cargo)

    def doProteins(self, fileiterator, datfile):
        """
        @brief reads the protein section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            isbound = self.isBoundary
            self.logs.datlog.info('Loading and Parsing proteins')

            while 1:
                allparams = {}
                while not isbound(self.cargo):
                    match = rx_protein.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')
                    self.cargo = fileiterator.next()
                datfile.addProteins(allparams)
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doProteins')
            raise

    def doQuery(self, fileiterator, datfile):
        """
        @brief reads the query sections and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        query = 'None'
        try:
            isbound = self.isBoundary
            while 1:
                allparams = {}
                match = rx_query.search(self.cargo)
                query = match.group('ID')
    #            if match.group('ID') == '1':
    #                self.logs.datlog.info('Loading and Parsing queries')
                allparams['query'] = query
                while not isbound(self.cargo):
                    match = rx_paramline.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')
                    self.cargo = fileiterator.next()

                datfile.addQuerySpectra(allparams)
                yield ('postquery', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'doQuery: last query = %s' % query)
            raise

    def doPostQuery(self, fileiterator, datfile):
        """
        @brief perfoms tasks after all the query data is loaded
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            datfile.doPostQuery()
            while 1:
                # Used to trigger the QC of the spectrum if required.
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doPostQuery')
            raise

    def doIndex(self, fileiterator, datfile):
        """
        @brief parses index section of the dat file
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing index')
            while 1:
                index = []
                while not isbound(self.cargo):
                    try:
                        key, value = self.cargo.split('=')
                        index.append((key, int(value)))
                    except:
                        pass
                    self.cargo = fileiterator.next()
                datfile.addIndex(index)
                yield ('boundary', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Called from doIndex')
            raise

    def doBoundary(self, fileiterator, datfile):
        """
        @brief reads fileiterator until it finds a boundry
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        while 1:
            while not rx_sectionh.search(self.cargo):
                self.cargo = fileiterator.next()
            match = rx_sectionh.search(self.cargo)
            sectionname = match.group('name')
            self.logs.datlog.debug("section detected: %s " % sectionname)
            if sectionname[:5] == 'query':
                sectionname = 'query'
            yield (sectionname, self.cargo)

    def doPostParsing(self, datfile):
        """
        @brief parses index section of the dat file
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        datfile.doPostParsing()


class DatParserPreCal(DatParser):

    def startParsing(self):
        """
        startParsing(string) -> integer

        Parses the dat file given by datfilepath.
        """
        datPath = str(self.datfileobj.filedic['datpath'])
        self.logs.datlog.info('Starting parsing file %s' % datPath)
        fin = open(datPath, 'r')

        try:
            self.cargo = fin.next()
        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'Empty dat file.')
            raise

        statedic = {'start': self.doSkip(fin, self.datfileobj),
                    'parameters': self.doSkip(fin, self.datfileobj),
                    'masses': self.doMasses(fin, self.datfileobj),
                    'unimod': self.doSkip(fin, self.datfileobj),
                    'enzyme': self.doSkip(fin, self.datfileobj),
                    'header': self.doSkip(fin, self.datfileobj),
                    'summary': self.doSummary(fin, self.datfileobj),
                    'et_summary': self.doSkip(fin, self.datfileobj),
                    'peptides': self.doPeptides(fin, self.datfileobj),
                    'postpeptides': self.doSkip(fin, self.datfileobj),
                    'et_peptides': self.doSkip(fin, self.datfileobj),
                    'proteins': self.doSkip(fin, self.datfileobj),
                    'query': self.doQuery(fin, self.datfileobj),
                    'postquery': self.doSkip(fin, self.datfileobj),
                    'index': self.doSkip(fin, self.datfileobj),
                    'boundary': self.doBoundary(fin, self.datfileobj)}

        coroutine = 'start'
        try:
            while 1:
                (coroutine, self.cargo) = statedic[coroutine].next()
        except StopIteration:
            self.doPostParsing(self.datfileobj)
            self.logs.datlog.info('Parsing complete')
        fin.close()

    def doPeptides(self, fileiterator, datfile):
        """
        @brief reads the peptide section and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        key = 'None'
        try:
            isbound = self.isBoundary

            self.logs.datlog.info('Loading and Parsing peptides')
            collection = {}
            lastquery = ''
            while 1:
                while not isbound(self.cargo):
                    match = rx_peptide.search(self.cargo)
                    if match:
                        value = match.group('value')
                        key = match.group('key')
                        if value == '-1':
                            # no matching peptide
                            datfile.stats['numspectra_nopeps'] += 1
                        else:
                            query = key.split('_')[0]
                            if query != lastquery:
                                lastquery = query
                                if collection:
                                    # add the query number to addpeptides
                                    datfile.addCalPeptide(collection)
                                    collection = {}
                            collection[key] = value
                    self.cargo = fileiterator.next()
                # add last peptide
                if collection:
                    # add the query number to addpeptides
                    datfile.addCalPeptide(collection)
                    pass
                yield ('postpeptides', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'doPeptides, last peptide = %s' % key)
            raise

    def doQuery(self, fileiterator, datfile):
        """
        @brief reads the query sections and puts the data into a dictionary, data is passed to the datfile object
        @param fileiterator <file object>: linked to the dat file
        @param datfile <datfile object>: containing all the dat file data and processing methods
        """
        query = 'None'
        try:
            isbound = self.isBoundary
            while 1:
                allparams = {}
                match = rx_query.search(self.cargo)
                query = match.group('ID')
    #            if match.group('ID') == '1':
    #                self.logs.datlog.info('Loading and Parsing queries')
                allparams['query'] = query
                while not isbound(self.cargo):
                    match = rx_paramline.search(self.cargo)
                    if match:
                        allparams[match.group('key')] = match.group('value')
                    self.cargo = fileiterator.next()

                datfile.addQuerySpectra(allparams)
                qry = int(query)
                datfile.spectra[qry]['spec_id'] = int(datfile.spectra[qry]['msmsid'][1:])
                datfile.spectra[qry]['rt'] = float(datfile.spectra[qry]['start'])
                yield ('postquery', self.cargo)

        except Exception, genEx:
            # catch exceptions and add some context data
            ExHa.addContext(genEx, 'doQuery: last query = %s' % query)
            raise
