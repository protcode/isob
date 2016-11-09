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

This file mainly deals with the handling of the command line data and reading
configuration files converting parameter values as directed.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob

"""

# python imports
from pathlib import Path
import sys
import numpy as np
import math
sys.path.insert(0, '..')

from CommonUtils.ConfigManager import ConfigManager
import CommonUtils.ExceptionHandler as ExHa
from CommonUtils.LoggingManager import Logger
from CommonUtils.hdf5Results import HDF5Results


class Top3quantController:
    def __init__(self,  hdf5quantprot, ncqcfg):

        self.hdf5File = hdf5quantprot.hdfFilePath
        self.hdf5quantprot = hdf5quantprot

        self.cfg = ncqcfg

    def dotop3quant(self):
        proteinset2top3 = {}
        top3peptoupdate = {}
        proteinset2pepdata = self.hdf5quantprot.getPeptideDataforTop3()
        proteinsalldata = self.hdf5quantprot.getSpecIDsUsedForQuant()
        ssdf = self.hdf5quantprot.getPeakRT_CS_lessthan30sfrompeak()
        for protein_group_no in proteinset2pepdata.keys():
            logs.log.debug('starting to process %s for top3' % protein_group_no)
            pcmdict = {}
            validsetids = [(x, y) for x, y in proteinset2pepdata[protein_group_no].iteritems()
                           if x in ssdf and x in proteinsalldata]
            if not validsetids:
                logs.log.debug('No valid quant data found for proteinset %s; will use all unique / fdr > threshold' %
                               protein_group_no)
                validsetids = [(x, y) for x, y in proteinset2pepdata[protein_group_no].iteritems()]
                logs.log.debug('there are %s validset ids from proteinset %s (not necessarily quantified)' %
                               (len(validsetids), protein_group_no))

            for spec_id, pepdata in validsetids:
                peptide = pepdata[0]
                logs.log.debug('processing spec_id %s with peptide data %s' % (spec_id, pepdata))
                try:
                    cstate = ssdf[spec_id][1]
                except KeyError:
                    logs.log.debug('no peak or charge_state data for spectrum %s' % spec_id)
                    continue

                pcmstring = '%s%s%s' % (peptide, cstate, pepdata[1])

                peak_intensity = ssdf[spec_id][0]

                rt = ssdf[spec_id][2]
                score = pepdata[2]
                logs.log.debug('pcmstring %s & intensity %s' % (pcmstring, peak_intensity))
                if pcmstring in pcmdict:
                    pcmdata = pcmdict[pcmstring]
                    if rt < pcmdata[0]:
                        logs.log.debug('start time less: use that one %s vs %s (score %s / %s)' %
                                       (rt, pcmdata[0], score, peptide))
                        pcmdict[pcmstring] = (rt, score, peak_intensity, peptide)
                    elif rt == pcmdata[0]:
                        if score > pcmdata[1]:
                            pcmdict[pcmstring] = (rt, score, peak_intensity, peptide)
                else:
                    pcmdict[pcmstring] = (rt, score, peak_intensity, peptide)
            peptide2peak_intensity = {}
            for rt, score, peak_intensity, peptide in pcmdict.values():
                logs.log.debug('processing pcm with RT %s, score %s, peak_intensity %s, peptide  %s' %
                               (rt, score, peak_intensity, peptide))
                if not peak_intensity:
                    continue
                try:
                    peptide2peak_intensity[peptide] += peak_intensity
                    logs.log.debug('adding for peptide %s (peak_intensity: %s) ' % (peptide, peak_intensity))
                except KeyError:
                    peptide2peak_intensity[peptide] = peak_intensity
                    logs.log.debug('first time for peptide %s (peak_intensity: %s) ' % (peptide, peak_intensity))
            top3 = [(y, x) for x, y in peptide2peak_intensity.iteritems()]
            top3.sort(key=lambda x: x[0], reverse=True)
            with np.errstate(invalid='ignore'):
                top3quantvalue = np.mean(np.array([math.log10(x[0]) for x in top3[:3]]))
            top3pepdata = [(x[1]+':'+str(math.log(x[0], 10))) for x in top3[:3]]
            top3peps = [x[1] for x in top3[:3]]
            proteinset2top3[protein_group_no] = top3quantvalue
            top3peptoupdate[protein_group_no] = top3peps
            logs.log.debug('top3 quant value %f %s (from %s )' %
                           (top3quantvalue, protein_group_no, ','.join(top3pepdata)))
        return proteinset2top3, top3peptoupdate


if __name__ == "__main__":
    logger = 0
    cfg = ConfigManager('./top3quant.cfg')
    performs2iupdate = 1
    try:
        configErr = cfg.evaluateCommandLineArgs(sys.argv)
        logParam = cfg.parameters['logging']
        logs = Logger(logParam['logfile'], logParam['loglevel'], logParam['screenlevel'])
        logs.setLog('top3quant')
        cfg.log = logs.log

        hdf5files = [Path(cfg.parameters['runtime']['datadir'], cfg.parameters['runtime']['filename'])]
        allpeptideratios = {}
        allsumionratiodata = {}
        counter = 0
        for hdf5Path in hdf5files:
            logs.log.info('hdf5Path %s' % hdf5Path)
            logs.log.info('starting')
            hdf5quanprot = HDF5Results(hdf5Path)
            hdf5quanprot.appendOpen()
            hdf5quanprot.cfg = cfg
            top3quant = Top3quantController(hdf5quanprot, cfg)
            top3datatoinsert, top3peptoupdate = top3quant.dotop3quant()
            hdf5quanprot.updateProteinHitTop3(top3datatoinsert)
            hdf5quanprot.updatePeptideTop3(top3peptoupdate)
            hdf5quanprot.addConfigParameters(cfg.parameters, 'postMascot', '3 top3quant')
            hdf5quanprot.close()

    except ExHa.czException as czEx:
        ExHa.reformatException(czEx)
        ExHa.addContext(czEx, 'Error during top3 run')
        ExHa.exportError2File(czEx, cfg.parameters['runtime']['datadir'] / Path('errors.error'))
        if logger:
            logger.log.warning(ExHa.oneLineRepr(czEx))
        else:
            print ExHa.multiLineRepr(czEx)

    except Exception as genEx:

        ExHa.reformatException(genEx)
        ExHa.addContext(genEx, 'Error during top3 run')
        ExHa.exportError2File(genEx, cfg.parameters['runtime']['datadir'] / 'errors.error')
        if logger:
            logger.log.warning(ExHa.oneLineRepr(genEx))
        else:
            print ExHa.multiLineRepr(genEx)
