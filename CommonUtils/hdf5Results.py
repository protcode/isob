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

This file handles the hdf5 interface for the results file data.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob

"""

from hdf5Base import hdf5Base
import sys

class HDF5Results:
    def __init__(self, hdfFilePath):
        self.hdfFilePath = hdfFilePath
        self.hdf = hdf5Base(hdfFilePath, forceCreation=True)

    def appendOpen(self):
        self.hdf.appendOpen()

    def readOpen(self):
        self.hdf.readOpen()

    def close(self):
        self.hdf.close()

    def createTables(self):
        tables = [('sample', 'ResultSample'), ('spectrum', 'ResultSpectrum'),
                  ('specquant', 'ResultSpecQuant'), ('peptide', 'ResultPeptide'),
                  ('proteinhit', 'ResultProteinHit'), ('statistics', 'ResultStatistics'),
                  ('proteinquant', 'ResultProteinQuant'), ('fdrdata', 'FDRData'),
                  ('configparameters', 'ResultConfigParameters')]
        for t in tables:
            self.hdf.createTable('/', t[0], t[1])

    def createQuantTables(self):
        """
        @brief : create tables to store quant data for each protein
        If already there then overwrite to prevent getting repeat data
        """
        tables = [('proteinquant', 'ResultProteinQuant')]
        for t in tables:
            if self.hdf.hasTable('/%s' % t[0]):
                self.hdf.createTable('/', t[0], t[1], delete=True)
        if not self.hdf.hasTable('/configparameters'):
            self.hdf.createTable('/', 'configparameters', 'ResultProcessingParameters')

    def createIndexes(self):
        """
        @brief organises the creation of indexes on appropriate tables in the hdf5 file
        @return:
        """

        hdf = self.hdf
        hdf.indexTable('/sample', ['sample_id'])
        hdf.indexTable('/spectrum', ['spectrum_id', 'sample_id'])
        hdf.indexTable('/specquant', ['spectrum_id', 'isotopelabel_id', 'protein_group_no'])
        hdf.completely_sorted_index_table('/peptide', 'protein_group_no')
        hdf.indexTable('/peptide', ['peptide_id', 'spectrum_id'])
        hdf.indexTable('/proteinhit', ['protein_group_no', 'protein_id'])
        hdf.indexTable('/proteinquant', ['protein_group_no'])

    def updateProteinHitTop3(self, top3data):
        hdf = self.hdf
        phitref = hdf.getTable('/proteinhit')
        for row in phitref:
            protein_group_no = row['protein_group_no']
            top3value = top3data[protein_group_no]
            row['top3'] = top3value
            row.update()
        phitref.flush()

    def updatePeptideTop3(self, top3peptoupdate):
        hdf = self.hdf
        pepref = hdf.getTable('/peptide')
        for row in pepref:
            protein_group_no = row['protein_group_no']
            top3peps = top3peptoupdate[protein_group_no]
            if row['peptide'] in top3peps:
                row['in_top3'] = 1
                row.update()
        pepref.flush()

    def writeSample(self, summary):
        write = []
        write.append(dict(sample_id=summary['sample_id'],
                          source_file=summary['source_file'].replace('.hdf5', '.raw'),
                          acquisition_time=summary['runtime'],
                          acquired_spectra=summary['acquired_spectra'],
                          mascot_matched_spectra=summary['mascot_matched_spectra'],
                          spectra_in_qc_proteins=summary['spectra_in_qc_proteins'],
                          quantified_spectra=summary['quantified_spectra'],
                          mean_precursor_ion_accuracy=summary['mean_precursor_ion_accuracy'],
                          sd_precursor_ion_accuracy=summary['sd_precursor_ion_accuracy'],
                          mean_reporter_ion_accuracy=summary['mean_reporter_ion_accuracy'],
                          sd_reporter_ion_accuracy=summary['sd_reporter_ion_accuracy'],
                          ))
        return self.hdf.appendRows('/sample', write)

    def writeSpectrum(self, spectra):
        write = []

        for data in spectra:
            write.append(dict(spectrum_id=data['spectrum_id'],
                              sample_id=data['sample_id'],
                              msms_id=int(data['msms_id'].replace('F', '')),
                              query=data['query'],
                              neutral_mass=data['neutral_mass'],
                              charge_state=data['charge_state'],
                              parent_ion=data['parent_ion'],
                              precursor_mz=data['precursor_mz'],
                              s2i=data['s2i'],
                              p2t=data['precursor_area'] / data['precursor_noise'],
                              rt=data['rt'],
                              survey_id=data['survey_id'],
                              start_time=data['start_time'],
                              peak_rt=data['peak_rt'],
                              peak_intensity=data['peak_inten_raw'],
                              peak_fwhm=data['peak_fwhm'],
                              sum_reporter_ions=data['sumreporterarea'],
                              quant_cancelled=0))
        return self.hdf.appendRows('/spectrum', write)

    def writeSpecQuant(self, specQuant):
        write = []
        for data in specQuant:
            write.append(dict(spectrum_id=data['spectrum_id'],
                              isotopelabel_id=data['isotopelabel_id'],
                              delta_seq=data['delta_seq'],
                              protein_group_no=data['protein_group_no'],
                              is_unique=data['is_unique'],
                              quant_raw=data['quant_raw'],
                              quant_isocorrected=data['quant_isocorrected'],
                              quant_allcorrected=data['quant_allcorrected'],
                              score=data['score'],
                              s2i=data['s2i'],
                              peptide_length=data['peptide_length'],
                              p2t=data['p2t'],
                              fdr_at_score=data['fdr_at_score'],
                              in_quantification_of_protein=data['in_quantification_of_protein'],
                              least_squares=data['least_squares'],
                              prior_ion_ratio=data['prior_ion_ratio'],
                              ms1source=data['ms1source']
                              ))
        return self.hdf.appendRows('/specquant', write)

    def writePeptide(self, peptides):
        write = []
        for data in peptides:
            write.append(dict(peptide_id=data['peptide_id'],
                              spectrum_id=data['spectrum_id'],
                              protein_group_no=data['protein_group_no'],
                              peptide=data['peptide'],
                              fdr_at_score=data['fdr_at_score'],
                              variable_modstring=data['variable_modstring'],
                              fixed_modstring=data['fixed_modstring'],
                              positional_modstring=data['positional_modstring'],
                              score=data['score'],
                              rank=data['rank'],
                              mw=data['mw'],
                              da_delta=data['da_delta'],
                              ppm_error=data['ppm_error'],
                              missed_cleavage_sites=data['missed_cleavage_sites'],
                              delta_seq=data['delta_seq'],
                              delta_mod=data['delta_mod'],
                              is_hook=data['is_hook'],
                              is_duplicate=data['is_duplicate'],
                              is_first_use_of_sequence=data['is_first_use_of_sequence'],
                              is_unique=data['is_unique'],
                              is_decoy=data['is_decoy'],
                              is_quantified=data['is_quantified'],
                              failed_fdr_filter=data['failed_fdr_filter'],
                              in_protein_inference=data['in_protein_inference'],
                              seq_start=data['seq_start'],
                              seq_end=data['seq_end']
                              ))
        return self.hdf.appendRows('/peptide', write)

    def writeFDRdata(self, fdrdata, data_type):

        write = []

        for score, data in fdrdata.iteritems():
            write.append(dict(score=score,
                              target_hits=data[0],
                              decoy_hits=data[1],
                              true_hits=data[2],
                              global_fdr=data[3],
                              local_fdr=data[4],
                              data_type=data_type))

        return self.hdf.appendRows('/fdrdata', write)

    def writeProteinHit(self, proteinHit):

        write = []
        for data in proteinHit:
            write.append(dict(protein_group_no=data['protein_group_no'],
                              protein_id=data['protein_id'],
                              gene_name=data['gene_name'],
                              description=data['description'],
                              mw=data['mw'],
                              ssm=data['ssm'],
                              total_score=data['total_score'],
                              hssm=data['hssm'],
                              hookscore=data['hookscore'],
                              upm=data['upm'],
                              is_decoy=data['is_decoy'],
                              protein_fdr=data['protein_fdr'],
                              max_score=data['max_score']
                              ))
        return self.hdf.appendRows('/proteinhit', write)

    def writeStatistics(self, statsDict):

        write = []
        keys = statsDict.keys()
        keys.sort()
        for s in keys:
            write.append(dict(statistic=s, value=statsDict[s]))
        return self.hdf.appendRows('/statistics', write)

    def writeProteinQuant(self, proteinQuantlist):

        return self.hdf.appendRows('/proteinquant', proteinQuantlist)

    def getFilterSpectra(self, whereclause):
        hdf = self.hdf
        tablePath = '/specquant'
        spectrumdata = hdf.getDataGeneral(tablePath, whereclause)
        self.cfg.log.debug('got %s records after %s' % (len(spectrumdata), whereclause))
        return set(spectrumdata['spectrum_id'])

    def getPeakRT_CS_lessthan30sfrompeak(self):
        hdf = self.hdf
        returndata = {}
        whereclause = 'start_time - peak_rt < 0.5'
        spectra = hdf.getDataGeneral('/spectrum', whereclause)
        for row in spectra:
            spectrum_id = row['spectrum_id']
            peak_inten = row['peak_intensity']
            charge_state = row['charge_state']
            peak_rt = row['peak_rt']
            returndata[spectrum_id] = (peak_inten, charge_state, peak_rt)

        return returndata

    def getSpecIDsUsedForQuant(self):
        hdf = self.hdf
        returndata = set()
        whereclause = 'in_quantification_of_protein==1'
        if hdf.hasTable('/specquant'):
            spectra = hdf.getDataGeneral('/specquant', whereclause)
            for row in spectra:
                spectrum_id = row['spectrum_id']
                returndata.add(spectrum_id)

        return returndata

    def getQuantData(self):
        hdf = self.hdf
        tablePath = '/specquant'
        ms2dataref = hdf.getTable(tablePath)
        return ms2dataref

    def updates2ivalues(self, s2icorrecteddata):
        """
        s2icorrecteddata
        :param s2icorrecteddata:is a dictionary of dictionaries, where keys are spectrum ids the second level is
          indexed by isotope label id
        :return:
        """

        hdf = self.hdf
        tablePath = '/specquant'
        ms2dataref = hdf.getTable(tablePath)
        for row in ms2dataref:
            correctedvalue = 0
            spectrum_id = row['spectrum_id']

            isotopelabel_id = row['isotopelabel_id']
            if spectrum_id in s2icorrecteddata:
                if isotopelabel_id in s2icorrecteddata[spectrum_id]:
                    correctedvalue = s2icorrecteddata[spectrum_id][isotopelabel_id]
                    self.cfg.log.debug('will update %s with corrected quantvalue %s' % (spectrum_id, correctedvalue))
            row['quant_allcorrected'] = correctedvalue
            row.update()
        ms2dataref.flush()

    def getProteinDatafromSpecQuant(self, proteinsgroupid):
        hdf = self.hdf
        tablePath = '/specquant'
        whereclause = 'protein_group_no == %s' % proteinsgroupid
        protdata = hdf.getDataGeneral(tablePath, whereclause)
        return protdata

    def getPeptidDataforSet(self, protein_group_no):
        hdf = self.hdf
        tablePath = '/peptide'
        whereclause = 'protein_group_no == %s' % protein_group_no
        pepdata = hdf.getDataGeneral(tablePath, whereclause)
        return pepdata

    def getPeptideDataforSets(self):
        hdf = self.hdf
        returndict = {}
        tablePath = '/peptide'
        pepdataref = hdf.getTable(tablePath)
        for row in pepdataref:
            protein_group_no = row['protein_group_no']
            peptide = row['peptide']
            spectrum_id = row['spectrum_id']
            try:
                returndict[protein_group_no][spectrum_id] = peptide
            except KeyError:
                returndict[protein_group_no] = {}
                returndict[protein_group_no][spectrum_id] = peptide
        return returndict

    def getAllProteinDatafromSpecQuant(self):
        hdf = self.hdf
        tablePath = '/specquant'
        whereclause = 'in_quantification_of_protein == 1'
        protdata = hdf.getDataGeneral(tablePath, whereclause)
        return protdata

    def updatein_quantification_of_protein(self, spectra, is_quantified_value):
        """
        updates the value of in_quantification_of_protein to value given by is_quantified_value
        ie for including or excluding a spectrum from quantification
        :param spectra:  list of spectrum ids to update
        :param is_quantified_value: value to which one should set the quantification of spectrum id

        :return:
        """

        hdf = self.hdf
        tablePath = '/specquant'
        specquantref = hdf.getTable(tablePath)
        refcounter = len(specquantref)
        for row in specquantref:
            spectrum_id = row['spectrum_id']
            if spectrum_id in spectra:
                row['in_quantification_of_protein'] = is_quantified_value
                row.update()
        specquantref.flush()

    def getQuantifiedPeptideCounts(self):
        hdf = self.hdf
        tablePath = '/peptide'
        whereclause = 'is_quantified== 1'
        potentially_quantfied_peptides = hdf.getDataGeneral(tablePath, whereclause)
        protein_group_no2quantpep = {}
        for data in potentially_quantfied_peptides:
            if not data['failed_fdr_filter']:
                protein_group_no = data['protein_group_no']
                try:
                    protein_group_no2quantpep[protein_group_no] += 1
                except KeyError:
                    protein_group_no2quantpep[protein_group_no] = 1
        return protein_group_no2quantpep

    def getProteinQuantData(self):
        hdf = self.hdf
        quantdata = hdf.readTable('/proteinquant')

        returndict = {}
        reference_label = 0

        for a in quantdata:
            protein_group_no = a['protein_group_no']
            isotopelabel_id = a['isotopelabel_id']
            reference_label = a['reference_label']
            fold_change = a['protein_fold_change']
            lower_conf = a['lower_confidence_level']
            upper_conf = a['upper_confidence_level']
            sum_quant_signal = a['sum_quant_signal']
            qssm = a['qssm']
            qupm = a['qupm']

            try:
                returndict[protein_group_no][isotopelabel_id] = dict(rel_fc=fold_change,
                                                                     delta_conf=upper_conf-lower_conf,
                                                                     signal_sum=sum_quant_signal,
                                                                     qupm=qupm,
                                                                     qssm=qssm)
            except KeyError:
                returndict[protein_group_no] = {}
                returndict[protein_group_no][isotopelabel_id] = dict(rel_fc=fold_change,
                                                                     delta_conf=upper_conf-lower_conf,
                                                                     signal_sum=sum_quant_signal,
                                                                     qupm=qupm,
                                                                     qssm=qssm)

        return returndict, reference_label

    def addConfigParameters(self, cfgParams, prePost, controller, sectionSkip=[]):

        if isinstance(sectionSkip, str):
            sectionSkip = [sectionSkip]
        elif not isinstance(sectionSkip, list):
            sectionSkip = list(sectionSkip)

        sectionSkip += ['runtime', 'logging']
        dataList = []
        for section in cfgParams:
            if section in sectionSkip:
                continue
            for param in cfgParams[section]:
                if param in ['paramconversion', '__name__']:
                    continue
                dataList.append(dict(analysis=prePost, application=controller, section=section, parameter=param,
                                     value=str(cfgParams[section][param])))
        if dataList:
            self.hdf.appendRows('/configparameters', dataList)

    def getPeptideDataforTop3(self):
        hdf = self.hdf
        returndict = {}
        tablePath = '/peptide'
        pepdataref = hdf.getTable(tablePath)
        for row in pepdataref:
            protein_group_no = row['protein_group_no']
            if protein_group_no not in returndict:
                returndict[protein_group_no] = {}
            peptide = row['peptide']
            positional_modstring = row['positional_modstring']
            rank = row['rank']
            score = row['score']
            spectrum_id = row['spectrum_id']
            is_unique = row['is_unique']
            failed_fdr_filter = row['failed_fdr_filter']

            if rank == 1 and is_unique and not failed_fdr_filter:

                try:
                    returndict[protein_group_no][spectrum_id] = (peptide, positional_modstring, score)
                except KeyError:
                    returndict[protein_group_no][spectrum_id] = (peptide, positional_modstring, score)

        return returndict

    def getH5DFQuantMeth(self):
        isotopes = self.hdf.readTable('/isotopes')
        if len(isotopes) == 0:
            return 0
        else:
            return isotopes[0]['method_id']



