"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2016 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive / Fusion
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This holds groups of peptides together with the protein accessions of proteins
they belong to. Stats about scores and counts are based on the peptide
members of these groups.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here:
https://github.com/protcode/isob
"""

# !/usr/bin/env python


# python imports
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, '..')

# CommonUtils imports
from CommonUtils.ConfigManager import ConfigManager
import CommonUtils.ExceptionHandler as ExHa
from CommonUtils.LoggingManager import Logger
from CommonUtils.hdf5Results import HDF5Results
import CommonUtils.progressbar as progBar
from CommonUtils.progressReport import progressReport


class CorrectS2IQuantController:
    def __init__(self,  nqchdf5corrects2iquant, ncqcfg):
        self.hdf5File = nqchdf5corrects2iquant.hdfFilePath
        self.hdf5corrects2iquant = nqchdf5corrects2iquant
        self.cfg = ncqcfg
        self.isotopecorrecteddata = {}
        self.s2icorrecteddata = {}
        self.fcratiodict = {}
        self.samplename = 'filepath'  # will be overwritten by real value.
        self.spectrumid2s2i = {}
        self.thisisotopecorrecteddata = {}

    def getQuantData(self):
        """@brief get ms2 data from hdf5 file, check whether needs to be isotope corrected or not (control) then
        pop it in a dictionary following isotope correction. Write back to HDF5 file if necessary
        """
        rawquanref = self.hdf5corrects2iquant.getQuantData()
        rawquan = {}
        spectrumid2s2i = {}

        for row in rawquanref.read():
            s2i = row['s2i']
            isotopelabel_id = row['isotopelabel_id']
            spectrum_id = row['spectrum_id']
            quant_isocorrected = row['quant_isocorrected']
            quant_raw = row['quant_raw']
            spectrumid2s2i[spectrum_id] = s2i
            if spectrum_id not in rawquan:
                rawquan[spectrum_id] = {}
            if quant_isocorrected:

                rawquan[spectrum_id][isotopelabel_id] = quant_isocorrected
            else:
                # use raw area
                rawquan[spectrum_id][isotopelabel_id] = quant_raw
        self.isotopecorrecteddata = rawquan
        self.spectrumid2s2i = spectrumid2s2i
        return rawquan

    def gets2icorrection(self, thisisotopecorrecteddata):
        """
        @brief get ratios one reporter ion / sum(reporter ion) for each
        """
        self.cfg.log.info('started gets2icorrection')
        ratiodict = {}
        sumdict = {}
        for areadata in thisisotopecorrecteddata.values():
            sumdata = sum(areadata.values())
            for isotoplabel_id, val in areadata.iteritems():
                if isotoplabel_id not in ratiodict:
                    ratiodict[isotoplabel_id] = []
                    sumdict[isotoplabel_id] = 0
                ratiodict[isotoplabel_id].append(val/float(sumdata))
                sumdict[isotoplabel_id] += val

        self.sumdict = sumdict
        self.cfg.log.debug('done gets2icorrection added sumdict %s ' % sumdict)
        return ratiodict

    def performS2Icorrection(self, correctionfactors):
        """
        @brief actully performs the s2i correction. The corrected data are written to the .hdf5 file in a separate
        column
        @correctionfactors dict normalized calculated factors giving amount of interfering signal from other reporter
        signals to be removed depending on s2i of ms1
        """
        spectrum2s2i = self.spectrumid2s2i
        self.cfg.log.info('we have %s s2i records' % len(spectrum2s2i))
        mys2icorrecteddata = {}

        progRep = progressReport(len(self.isotopecorrecteddata), self.hdf5File.stem, 'S2I correcting', 'precursors')

        for spectrum_id, data in self.isotopecorrecteddata.iteritems():
            self.cfg.log.debug('spectrum id is %s' % spectrum_id)
            # perform correction only if there is actually an S2I value for that spectrum.
            if spectrum_id in spectrum2s2i:
                progRep.next()
                s2ivalue = round(spectrum2s2i[spectrum_id], 3)
                self.cfg.log.debug('spectrum_id for s2i correction %s, s2i value %s ' % (spectrum_id, s2ivalue))

                for isotopelabel_id, val in data.iteritems():

                    s2icorr = val - ((1 - s2ivalue) * correctionfactors[isotopelabel_id] * sum(data.values()))
                    self.cfg.log.debug('calculated s2i corrected value %f from original %f for %s ' %
                                       (s2icorr, val, isotopelabel_id))
                    if s2icorr < 0:  # no signal can be less than zero
                        s2icorr = 0
                    if spectrum_id not in mys2icorrecteddata:
                        mys2icorrecteddata[spectrum_id] = {}
                    mys2icorrecteddata[spectrum_id][isotopelabel_id] = s2icorr
            else:
                self.cfg.log.debug('no spectrum_id (%s) for s2i correction ' % spectrum_id)
                for isotopelabel_id, val in data.iteritems():
                    if spectrum_id not in mys2icorrecteddata:
                        mys2icorrecteddata[spectrum_id] = {}
                    mys2icorrecteddata[spectrum_id][isotopelabel_id] = val

        progRep.endReport()
        self.s2icorrecteddata = mys2icorrecteddata
        self.cfg.log.info('done performS2Icorrection')
        return mys2icorrecteddata


if __name__ == "__main__":
    logger = 0
    cfg = ConfigManager('./corrects2iquant.cfg')

    try:
        configErr = cfg.evaluateCommandLineArgs(sys.argv)
        dataDir = cfg.parameters['runtime']['datadir']
        # start the logging process
        logParam = cfg.parameters['logging']
        logPath = Path(dataDir.joinpath(logParam['logdir']))
        if not logPath.exists():
            logPath.mkdir(parents=True)
        logFile = logPath.joinpath(logParam['logfile'])
        logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], False)
        logger.setLog('corrects2iquant')
        cfg.log = logger.log
        hdf5files = [Path(cfg.parameters['runtime']['datadir'], cfg.parameters['runtime']['filename'])]
        allpeptideratios = {}
        allsumionratiodata = {}
        corrects2iquantoblist = []
        counter = 0
        for hdf5Path in hdf5files:
            logger.log.info('hdf5Path %s' % hdf5Path)
            logger.log.info('starting')
            hdf5corrects2iquant = HDF5Results(hdf5Path)
            hdf5corrects2iquant.appendOpen()
            hdf5corrects2iquant.samplename = hdf5Path.name
            hdf5corrects2iquant.cfg = cfg
            hdf5corrects2iquant.addConfigParameters(cfg.parameters, 'postMascot', '4 corrects2iquant')

            s2iq = CorrectS2IQuantController(hdf5corrects2iquant, cfg)
            counter += 1
            s2iq.samplename = '0%s' % counter
            corrects2iquantoblist.append(s2iq)
            isotopecorrecteddata = s2iq.getQuantData()

            sumionratiodict = s2iq.gets2icorrection(isotopecorrecteddata)
            for isotopelabel_id, data in sumionratiodict.iteritems():
                if isotopelabel_id not in allsumionratiodata:
                    allsumionratiodata[isotopelabel_id] = []
                allsumionratiodata[isotopelabel_id].extend(data)

        allfractionbgratios = {}
        for isotopelabel_id, data in allsumionratiodata.iteritems():
            allfractionbgratios[isotopelabel_id] = np.median(data)
        # second normalization so that the bg-ratios all add to 1
        starting_allfractionbgratios = allfractionbgratios.copy()
        for isotopelabel_id, data in allfractionbgratios.iteritems():
            allfractionbgratios[isotopelabel_id] = data / sum(starting_allfractionbgratios.values())
        logger.log.debug(('allfractionbgratios are %s' % str(allfractionbgratios)))
        for corrects2iquantob in corrects2iquantoblist:
            # perform correction for each of the analyzed .hdf5 files.
            s2icorrecteddata = corrects2iquantob.performS2Icorrection(allfractionbgratios)
            corrects2iquantob.hdf5corrects2iquant.updates2ivalues(s2icorrecteddata)
            hdf5corrects2iquant.close()

    except ExHa.czException as czEx:
        ExHa.reformatException(czEx)
        ExHa.addContext(czEx, 'Error during corrects2iquant run')
        ExHa.exportError2File(czEx, cfg.parameters['runtime']['datadir'] / Path('errors.error'))
        if logger:
            logger.log.warning(ExHa.oneLineRepr(czEx))
        else:
            print ExHa.multiLineRepr(czEx)

    except Exception as genEx:

        ExHa.reformatException(genEx)
        ExHa.addContext(genEx, 'Error during corrects2iquant run')
        ExHa.exportError2File(genEx, cfg.parameters['runtime']['datadir'] / 'errors.error')
        if logger:
            logger.log.warning(ExHa.oneLineRepr(genEx))
        else:
            print ExHa.multiLineRepr(genEx)
