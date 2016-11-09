"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2016 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive / Fusion
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This module takes data from a results.hdf5 file and prints the data to one
of 3 text files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here:
https://github.com/protcode/isob
"""


import os
import sys
from pathlib import Path
import numpy as np
np.seterr(all='ignore')  # seterr to ignore
import re

RE_MOD = re.compile('(.+):(\d+)$')

from CommonUtils.hdf5Results import HDF5Results
from CommonUtils.ConfigManager import ConfigManager
from CommonUtils.LoggingManager import Logger
from CommonUtils.tools import Misc
import CommonUtils.progressbar as progBar
import CommonUtils.QuantMethodHandler as quantHandler
from CommonUtils.progressReport import progressReport


def getUniprotData():
    id2updata = {}
    try:
        fh = open('essentialuniprotdata.txt')
        for line in fh:
            protein_id, desc, gene = line.strip().split('\t')
            id2updata[protein_id] = (desc, gene)
    except:
        logger.log.info('issue with essentialuniprotdata.txt file: some info may be missing in the output')
    return id2updata


def outputResults(hdf5file):
    logger.log.info('processing file: %s' % hdf5file.name)
    hdf = HDF5Results(hdf5file)
    hdf.readOpen()
    sample2source = collectSummaryData(hdf)
    collectPeptideData(hdf, sample2source)
    collectProteinData(hdf)

    hdf.close()


def renameFile(file, suffix):
    newFile = file.parent.joinpath(file.stem + suffix + '.txt')
    return newFile


def makemodificationsstring(positional_modstring, peptide):
    """
    @brief: simple method to add extra residue information to positional modstring
    ie Carbamidomethyl:4; TMT6plex:5 --> C4 Carbamidomethyl; K5 TMT6plex
    GVACK
    positional_modstring:
    """
    modificdations = ''
    template = '%s %s%s; '
    template_term = '%s %s; '
    pmods = positional_modstring.split(';')
    for pmod in pmods:
        try:
            if RE_MOD.search(pmod):
                mod = RE_MOD.search(pmod).group(1)
                pos = RE_MOD.search(pmod).group(2)
            peptidepos = int(pos) - 1
            if peptidepos < 0:
                residue = 'N-term'
                modificdations += template_term % (mod, residue)

            elif peptidepos > len(peptide):
                residue = 'C-term'
                modificdations += template_term % (mod, residue)
            else:
                residue = peptide[peptidepos]
                modificdations += template % (mod, residue, pos)
        except ValueError:
            print 'error processing positional modstring %s' % positional_modstring
            return 'n/d'
    return modificdations


def preparePeptideData(pep, proteins, fdr_data, specs, quant):
    protein = proteins[pep['protein_group_no']]
    specID = pep['spectrum_id']
    fdr_at_score = fdr_data[round(pep['score'])]
    try:
        in_top3 = pep['in_top3']
    except KeyError:
        in_top3 = 'NA'
    delta_mod = round(pep['delta_mod']) if pep['delta_mod'] > -1 else 'NA'

    peptideData = dict(sequence=pep['peptide'],
                       rank=pep['rank'],
                       mw=pep['mw'],
                       ppm_error=pep['ppm_error'],
                       is_unique=pep['is_unique'], score=pep['score'],
                       is_decoy=pep['is_decoy'],
                       failed_fdr_filter=pep['failed_fdr_filter'],
                       protein_group_no=pep['protein_group_no'],
                       protein_ids=protein,
                       in_protein_inference=pep['in_protein_inference'],
                       fdr_at_score=fdr_at_score,
                       delta_mod=delta_mod,
                       in_top3=in_top3,
                       seq_start=pep['seq_start'],
                       seq_end=pep['seq_end'],
                       positional_modstring=pep['positional_modstring'],
                       prior_ion_ratio='NA',
                       least_squares='NA',
                       ms1source='NA'
                       )

    peptideData.update(specs[specID])

    if peptideData['positional_modstring']:
        peptideData['modifications'] = makemodificationsstring(peptideData['positional_modstring'],
                                                               peptideData['sequence'])
    else:
        peptideData['modifications'] = ''
    if specID in quant:
        peptideData.update(quant[specID])
        if {'prior_ion_ratio', 'least_squares', 'ms1source'} & set(quant[specID].keys()):
            peptideData['prior_ion_ratio'] = quant[specID]['prior_ion_ratio']
            peptideData['least_squares'] = quant[specID]['least_squares']
            peptideData['ms1source'] = quant[specID]['ms1source']
    else:
        peptideData['in_quantification_of_protein'] = 0

    peptideData['proteins'] = '|'.join(peptideData['protein_ids'])
    return peptideData


def collectPeptideData(hdfObject, sample2source):
    logger.log.info('generating Peptide based output')
    hdf = hdfObject.hdf
    outputFile = renameFile(hdf.filePath, '_peptides')

    # extract required protein data
    logger.log.info('loading protein data')
    proteinhit = hdf.readTable('/proteinhit')
    proteins = {}
    for prot in proteinhit:
        try:
            proteins[prot['protein_group_no']].append(prot['protein_id'])
            proteins[prot['protein_group_no']].sort()
        except KeyError:
            proteins[prot['protein_group_no']] = [prot['protein_id']]
    proteinhit = None

    # extract required spectrum data
    logger.log.info('loading spectrum data')
    spectra = hdf.readTable('/spectrum')
    specs = {}
    for sp in spectra:
        source_file = sample2source[sp['sample_id']]
        specs[sp['spectrum_id']] = dict(source_file=source_file,
                                        msms_id=sp['msms_id'],
                                        charge_state=sp['charge_state'],
                                        precursor_mz=sp['precursor_mz'],
                                        retention_time = sp['rt'],
                                        peak_intensity=sp['peak_intensity'],
                                        s2i=sp['s2i'],
                                        p2t=sp['p2t'])

    # extract quantification data
    logger.log.info('loading quantification data')
    specquant = hdf.readTable('/specquant')
    quant = {}
    usedIsotopes = set()
    addMS1data = False
    if set(cfg.parameters['general']['ms1params']) & set([x for x in specquant.dtype.names]):
        addMS1data = True
    for sq in specquant:
        id = sq['spectrum_id']
        if id not in quant:
            quant[id] = dict(in_quantification_of_protein=sq['in_quantification_of_protein'])
        quant[id][sq['isotopelabel_id']] = sq['quant_allcorrected']
        if addMS1data:
            quant[id]['prior_ion_ratio'] = sq['prior_ion_ratio']
            quant[id]['least_squares'] = sq['least_squares']
            quant[id]['ms1source'] = sq['ms1source']
        usedIsotopes.add(sq['isotopelabel_id'])

    usedIsotopes = sorted(usedIsotopes)
    outString = 'protein_group_no\tprotein_id\tsequence\tmodifications\tmw'
    outString += '\tprecursor_mz\tcharge_state\tppm_error\tscore\tfdr_at_score\tdelta mod\trank\tmsms_id\tsource_file'
    outString += '\tretention time\tpeak_intensity\ts2i\tp2t\tis_unique\tin_quantification_of_protein\tin_protein_inference'
    outString += '\tin_top3\tseq_start\tseq_end\tis_decoy'
    if addMS1data:
        outString += '\tprior ion ratio\tleast squares\tms1 source'
    if usedIsotopes:
        isotope_data = dict([(str(i), i) for i in usedIsotopes])
        try:
            y = quantHandler.QuantMethods()
            g = y.getMethodByIsotope(usedIsotopes[0])
            for id, data in g['quantmasses'].iteritems():
                isotope_data[id] = data[0]['name']
        except:
            print 'error getting label name data, just using names present in .hdf5 file'
        for iso in usedIsotopes:
            outString += '\tsig_%s' % isotope_data[iso]
    # open text file output
    f_out = open(str(outputFile), 'w')
    f_out.write(outString + '\n')
    # integrate other data with peptide data and output to text file
    logger.log.info('loading peptide data')
    out_string_template = '%(protein_group_no)i\t%(proteins)s\t%(sequence)s\t%(modifications)s\t%(mw)f\t'
    out_string_template += '%(precursor_mz)f\t%(charge_state)i\t%(ppm_error)f\t%(score).0f\t%(fdr_at_score).3f\t'
    out_string_template += '%(delta_mod)s\t%(rank)i\t%(msms_id)i\t%(source_file)s\t%(retention_time)f\t'
    out_string_template += '%(peak_intensity)f\t%(s2i)f\t%(p2t)f\t'
    out_string_template += '%(is_unique)i\t%(in_quantification_of_protein)i\t'
    out_string_template += '%(in_protein_inference)f\t%(in_top3)s\t%(seq_start)s\t%(seq_end)s\t%(is_decoy)s'
    if addMS1data:
        out_string_template += '\t%(prior_ion_ratio)s\t%(least_squares)s\t%(ms1source)s'
    fdr_data = dict([(x['score'], x['global_fdr']) for x in hdf.readTable('/fdrdata') if x['data_type'] == 'peptide'])
    peptidetable = hdf.getTable('/peptide')  # get reference to peptide table on disk
    current_pepgroupid = None
    tmplist = []

    progRep = progressReport(len(peptidetable), hdf.filePath.stem, 'load peptides', 'peptides')
    for idx, p in enumerate(peptidetable.itersorted('protein_group_no')):
        # pBar.update(idx)
        progRep.report(idx)
        if current_pepgroupid == p['protein_group_no']:
            tmplist.append(preparePeptideData(p, proteins, fdr_data, specs, quant))
        else:
            if tmplist:
                tmplist = sorted(tmplist, key=lambda y: y['seq_start'])
                for x in tmplist:
                    outString = out_string_template % x
                    for iso in usedIsotopes:
                        try:
                            outString += '\t%f' % x[iso]
                        except KeyError:
                            outString += '\tNA'

                    f_out.write(outString + '\n')
            tmplist = [preparePeptideData(p, proteins, fdr_data, specs, quant)]
        current_pepgroupid = p['protein_group_no']

    progRep.endReport()
    # not to forget the last protein groups data in tmplist
    if tmplist:
        tmplist = sorted(tmplist, key=lambda y: y['seq_start'])
        for x in tmplist:
            outString = out_string_template % x
            for iso in usedIsotopes:
                try:
                    outString += '\t%f' % x[iso]
                except KeyError:
                    outString += '\tNA'

            f_out.write(outString + '\n')
    f_out.close()
    return


def collectProteinData(hdfObject):
    logger.log.info('processing Protein data')
    # extract required protein data
    hdf = hdfObject.hdf
    outputFile = renameFile(hdf.filePath, '_proteins')
    logger.log.info('loading protein data')
    proteinhit = hdf.readTable('/proteinhit')
    allquantdata, referenceid = hdfObject.getProteinQuantData()
    potentially_quant_peps = hdfObject.getQuantifiedPeptideCounts()
    proteins = {}
    try:
        usedIsotopes = sorted(allquantdata.values()[0].keys())
    except:
        usedIsotopes = []
    translatedData = dict()
    translatedData[0] = ''
    referenceid = int(referenceid)
    essentialUPdata = getUniprotData()

    if usedIsotopes:
        translatedData = dict([(i, i) for i in usedIsotopes])
        try:
            y = quantHandler.QuantMethods()
            g = y.getMethodByIsotope(usedIsotopes[0])
            for id, data in g['quantmasses'].iteritems():
                translatedData[id] = data[0]['name']
        except:
            print 'error getting label name data, just use .hdf5 file content'
    for prot in proteinhit:
        protein_id = prot['protein_id']
        ssm = prot['ssm']
        protein_group_no = prot['protein_group_no']
        try:
            pqssm = potentially_quant_peps[protein_group_no]
        except KeyError:
            pqssm = 0

        is_decoy = prot['is_decoy']
        protein_fdr = prot['protein_fdr']
        try:
            top3 = prot['top3']
        except IndexError:
            top3 = np.nan
        if np.isnan(top3):
            top3 = 'NA'
        try:
            quantdata = allquantdata[protein_group_no]
            qupm = max([x['qupm'] for x in quantdata.values()])
            qssm = max([x['qssm'] for x in quantdata.values()])
        except KeyError:
            quantdata = []
            qupm = 0
            qssm = 0

        description = prot['description']
        gene_name = prot['gene_name']
        if gene_name == 'None' or description == 'missing':
            try:
                if protein_id in essentialUPdata:
                    description, gene_name = essentialUPdata[protein_id]
                else:
                    # can try splitting on the - in case the uniprot did not include isoforms
                    splitprotein_id = protein_id.split('-')[0]
                    description, gene_name = essentialUPdata[splitprotein_id]
            except KeyError:
                pass
        mw = prot['mw']
        total_score = prot['total_score']
        upm = prot['upm']
        max_score = prot['max_score']

        try:
            proteins[protein_group_no]['protein_id'].append(protein_id)
            proteins[protein_group_no]['description'].append(description)
            proteins[protein_group_no]['gene_name'].append(gene_name)
            proteins[protein_group_no]['mw'].append(mw)
            proteins[protein_group_no]['total_score'] = total_score
            proteins[protein_group_no]['quantdata'] = quantdata
            proteins[protein_group_no]['ssm'] = ssm
            proteins[protein_group_no]['qupm'] = qupm
            proteins[protein_group_no]['qssm'] = qssm
            proteins[protein_group_no]['upm'] = upm
            proteins[protein_group_no]['fqssm'] = pqssm - qssm
            proteins[protein_group_no]['is_decoy'].add(is_decoy)
            proteins[protein_group_no]['protein_fdr'] = protein_fdr
            proteins[protein_group_no]['max_score'] = max_score
            proteins[protein_group_no]['top3'] = top3
        except KeyError:
            proteins[protein_group_no] = {}
            proteins[protein_group_no]['protein_id'] = [protein_id]
            proteins[protein_group_no]['quantdata'] = quantdata
            proteins[protein_group_no]['description'] = [description]
            proteins[protein_group_no]['gene_name'] = [gene_name]
            proteins[protein_group_no]['mw'] = [mw]
            proteins[protein_group_no]['total_score'] = total_score
            proteins[protein_group_no]['ssm'] = ssm
            proteins[protein_group_no]['upm'] = upm
            proteins[protein_group_no]['qssm'] = qssm
            proteins[protein_group_no]['qupm'] = qupm
            proteins[protein_group_no]['fqssm'] = pqssm - qssm
            proteins[protein_group_no]['is_decoy'] = {is_decoy}
            proteins[protein_group_no]['protein_fdr'] = protein_fdr
            proteins[protein_group_no]['max_score'] = max_score
            proteins[protein_group_no]['top3'] = top3
    f_out = open(str(outputFile), 'w')
    outString = 'protein_group_no\tprotein_id\tdescription\tgene_name\tmw\ttop3\tprotein_fdr\tmax_score\ttotal_score\t'
    outString += 'ssm\tupm\tfqssm\tqssm\tqupm\n'
    quantfieldnames = ['signal_sum', 'rel_fc', 'delta_conf']
    if usedIsotopes:
        outString = outString[:-1] + '\treference_label\t'
        for columnname in quantfieldnames:
            for iso in usedIsotopes:
                title = columnname + '_%s' % translatedData[iso]
                outString += title+'\t'
        outString = outString[:-1] + '\n'

    f_out.write(outString)
    ignore_reverse_hits = cfg.parameters['general']['ignore_reverse_hits']
    MiscTools = Misc()
    for protein_group_no in proteins:
        proteindata = proteins[protein_group_no]
        if len(proteindata['protein_id']) > 1:
            sortedproteinacs = sorted(proteindata['protein_id'])
            sortedgenedata = []
            sortedmwdata = []
            sorteddescdata = []
            new2oldorder = {}
            for idx, spac in enumerate(sortedproteinacs):
                oldidx = proteindata['protein_id'].index(spac)
                new2oldorder[idx] = oldidx
            for new, old, in new2oldorder.iteritems():

                sortedgenedata.append(proteindata['gene_name'][old])
                sortedmwdata.append(proteindata['mw'][old])
                sorteddescdata.append(proteindata['description'][old])
            proteindata['protein_id'] = sortedproteinacs
            proteindata['gene_name'] = MiscTools.uniqifyList(sortedgenedata)
            proteindata['description'] = sorteddescdata
            proteindata['mw'] = sortedmwdata

        if ignore_reverse_hits and min(proteindata['is_decoy']) == 1:
            # do not print out the reverse hits.
            continue

        if len(set(proteindata['gene_name'])) == 1:
            #  can convert to a set if only 1 gene present, otherwise keep ALL entries.
            proteindata['gene_name'] = set(proteindata['gene_name'])

        proteindata2out = {}
        for k, v in proteindata.iteritems():
            if isinstance(v, list) or isinstance(v, set):
                proteindata2out[k] = '|'.join([str(x).replace('None', 'NA') for x in v])
            else:
                proteindata2out[k] = v
        proteindata2out['referenceid'] = translatedData[referenceid]
        quantdata = proteindata['quantdata']
        proteindata2out['protein_group_no'] = protein_group_no
        outString = '%(protein_group_no)i\t%(protein_id)s\t%(description)s\t%(gene_name)s\t%(mw)s\t%(top3)s\t'
        outString += '%(protein_fdr)s\t%(max_score).0f\t%(total_score).0f\t%(ssm)s\t'
        outString += '%(upm)s\t%(fqssm)s\t%(qssm)s\t%(qupm)s\t%(referenceid)s'
        outString = outString % proteindata2out
        for idx in quantfieldnames:
            for iso in usedIsotopes:

                try:
                    data = str(quantdata[iso][idx]).replace('-1.0', 'NA')
                    outString += '\t%s' % data
                except KeyError:
                    outString += '\tNA'
                except IndexError:
                    outString += '\tNA'
        f_out.write(outString + '\n')

    f_out.close()


def collectSummaryData(hdfObject):
    hdf = hdfObject.hdf
    sample2source = {}
    outputFile = renameFile(hdf.filePath, '_summary')
    f_out = open(str(outputFile), 'w')

    titles = ('sample_id', 'source_file', 'acquired_spectra', 'mascot_matched_spectra', 'spectra_in_qc_proteins',
              'quantified_spectra', 'mean_precursor_ion_accuracy', 'sd_precursor_ion_accuracy',
              'mean_reporter_ion_accuracy', 'sd_reporter_ion_accuracy')

    f_out.write('%s\n' % '\t'.join(titles))

    # load and output the individual sample data
    sampleSummary = hdf.readTable('/sample')
    sampleSummary = sorted(sampleSummary, key=lambda x: x['sample_id'])

    for ss in sampleSummary:
        outList = []
        for ti in titles:
            outList.append(str(ss[ti]))
        sample2source[ss['sample_id']] = ss['source_file']

        f_out.write('%s\n' % '\t'.join(outList))

    # load and output the combined summary data
    statistics = hdf.readTable('/statistics')

    combinedData = {}
    for s in statistics:
        combinedData[s['statistic']] = s['value']

    combinedData['sample_id'] = 'all'
    combinedData['source_file'] = 'combined'
    outList = []
    for ti in titles:
        outList.append(combinedData[ti])

    f_out.write('%s\n' % '\t'.join(outList))

    # create separate table for the protein numbers
    countTarget = 0
    countDecoy = 0
    for i in range(6):
        countTarget += int(combinedData['%i hooks target_protein_hits' % i])
        countDecoy += int(combinedData['%i hooks decoy_protein_hits' % i])

    f_out.write('\n\ntype\tprotein_hits\n')
    f_out.write('target\t%i\n' % countTarget)
    f_out.write('decoy\t%i\n' % countDecoy)

    # create a new table for peptide FDR levels
    fdr = hdf.readTable('/fdrdata')
    fdr.sort(order=['data_type', 'score'])

    f_out.write('\n\nfdr\tmin_mascot_score\n')
    fdrThresholds = {}
    for thresh in [5.0, 1.0]:
        filter = (fdr['data_type'] == 'peptide') & (fdr['global_fdr'] <= thresh / 100)
        fdrFilt = fdr[filter]
        fdrThresholds[thresh] = min(fdrFilt['score'])
        f_out.write('%i%%\t%i\n' % (thresh, fdrThresholds[thresh]))

    # create table for parameters
    param = hdf.readTable('/configparameters')
    param = sorted(param, key=lambda x: (x['application'], x['section']))
    param = sorted(param, key=lambda x: x['analysis'], reverse=True)

    f_out.write('\n\nanalysis\tapplication\tsection\tparameter\tvalue\n')
    for pa in param:
        f_out.write('%(analysis)s\t%(application)s\t%(section)s\t%(parameter)s\t%(value)s\n' % pa)

    f_out.close()
    return sample2source


if __name__ == '__main__':

    configPath = './outputResults.cfg'
    cfg = ConfigManager(configPath)
    ret = cfg.evaluateCommandLineArgs(sys.argv)

    pID = os.getpid()
    dataDir = cfg.parameters['runtime']['datadir']
    logParam = cfg.parameters['logging']
    logPath = Path(dataDir.joinpath(logParam['logdir']))
    if not logPath.exists():
        logPath.mkdir(parents=True)
    logFile = logPath.joinpath(logParam['logfile'])

    logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], False)
    logger.setLog('outputResults')

    logger.log.info('started outputResults application')
    filefiltervalue = cfg.parameters['runtime']['filefilter']
    # this weird construct makes sure that we only try and get outputs from .hdf5 files
    # if hdf5 alrady in filefilter name then it's removed and readded, otherwise it's added
    filefiltervalue = filefiltervalue.replace('.hdf5', '') + '.hdf5'
    for hdf in dataDir.glob(filefiltervalue):

        outputResults(hdf)

    logger.log.info('finished outputResults')
