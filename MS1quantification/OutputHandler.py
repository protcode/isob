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

# CommonUtils imports
import CommonUtils.ExceptionHandler as ExHa

# project imports


class OutputHandler:
    """
    @brief class to read and write pyMSsafe data into HDF5 format
    """
    def __init__(self, hdf5file, cfg, importGroup):
        self.hdf5 = hdf5file
        self.cfg = cfg
        self.importGroup = importGroup

        self.useSecondary = cfg.parameters['models']['compare']

        if self.cfg.parameters['models']['primary_model'].lower() == 'exact':
            self.primary = 'exact'
            self.secondary = 'averagine'
        else:
            self.primary = 'averagine'
            self.secondary = 'exact'


    def writeIsotopeData(self, specid, theory, label2id, isotopes):
        """
        appends the isotope data to the quanIso table
        @param specid <integer>: the spectrum identifier for the MS/MS event leading to the peptide identification
        @param theory <dictionary>: containing the theoretical isotope data (from the elemental composition)
        @param label2id  <dictionary>: containing the mapping of the label names to the isotope IDs
        @param isotopes <dictionary>: containing the measured isotope data
        @return:
        """
        try:
            tabGroup = self.importGroup
            # tabGroup = 'rawdata'
            tabName = 'quaniso'
            tableLoc = '/%s/%s' % (tabGroup, tabName)
            tabDef = 'QuanIso'

            rowsToAdd = []

            if not self.hdf5.hasTable(tableLoc):
                ok = self.hdf5.createTable(tabGroup, tabName, tabDef)

            colsUsed = dict(spec_id=1, isolabel_id=1, source=1, isointen_m1=1, isointen_0=1, isointen_p1=1)

            for meth in isotopes:
                for label in theory['byIncreasingMZ']:

                    trow = {}

                    if isotopes[meth][label]:
                        used = isotopes[meth][label]['cluster']['usedInten']

                        trow['spec_id'] = specid
                        trow['isolabel_id'] = label2id[label]
                        trow['source'] = meth
                        if -1 in used:
                            trow['isointen_m1'] = used[-1]
                        trow['isointen_0'] = used[0]
                        trow['isointen_p1'] = used[1]

                        for i in range(2, 6):
                            if i in used:
                                k = 'isointen_p' + str(i)
                                trow[k] = used[i]
                                colsUsed[k] = 1

                        rowsToAdd.append(trow)

            for label in theory['byIncreasingMZ']:
                trow = {}
                used = theory[label]['primary']
                trow['spec_id'] = specid
                trow['isolabel_id'] = label2id[label]
                trow['source'] = 'theory'
                if -1 in used:
                    trow['isointen_m1'] = used[-1]
                trow['isointen_0'] = used[0]
                trow['isointen_p1'] = used[1]

                for i in range(2, 6):
                    if i in used:
                        k = 'isointen_p' + str(i)
                        trow[k] = used[i]
                        colsUsed[k] = 1

                rowsToAdd.append(trow)

            ok = self.hdf5.appendRows(tableLoc, rowsToAdd, colsUsed.keys())
        except Exception as genEx:
            ExHa.addContext(genEx, 'Error during OutputHandler.writeIsotopeData')
            raise

    def writeQuant(self, quant):
        """
        @brief writes quantification data to the quant table and additional MS1 quant data to the quanextra table
        @param quant <dictionary>: contains all the quantification data
        @return:
        """

        # create quan table if needed
        tabGroup = self.importGroup
        # tabGroup = 'rawdata'
        tabName = 'quan'
        tableLoc = '/%s/%s' % (tabGroup, tabName)
        tableDef = 'Quan'

        # create main quan data table
        if not self.hdf5.hasTable(tableLoc):
            ok = self.hdf5.createTable(tabGroup, tabName, tableDef)

        # add data to the quan table
        quanRow = {}
        quanCols = ['spec_id', 'isolabel_id', 'area', 'corrected', 'inten', 'survey_id', 'mzdiff', 'ppm']

        for col in quanCols:
            quanRow[col] = quant[col]

        ok = self.hdf5.appendRows(tableLoc, [quanRow], quanCols)

        # create quanextra table
        tabName = 'quanextra'
        tableDef = 'QuanExtra'
        tableLoc = '/%s/%s' % (tabGroup, tabName)

        if not self.hdf5.hasTable(tableLoc):
            ok = self.hdf5.createTable(tabGroup, tabName, tableDef)


        quanExtraRow = dict(spec_id=quant['spec_id'],
                            ms1SpecID=quant['ms1SpecID'],
                            isolabel_id=quant['isolabel_id'],
                            theor_mz=quant['monoMZ'],
                            sum_theor_int=quant['sum_theor_int'],
                            source=quant['source'],
                            prior_ion_ratio_xic=quant['-1ratio_xic'],
                            prior_ion_ratio_apex=quant['-1ratio_apex'],
                            prior_ion_ratio_prior=quant['-1ratio_prior'],
                            prior_ion_ratio_final=quant['-1ratio_final'],
                            sum_least_squares=quant['cl_sum_fit'],
                            max_least_squares=quant['cl_max_fit'],
                            xic_mz=quant['XIC']['mz'],
                            xic_rt=quant['XIC']['rt'],
                            xic_least_squares=quant['XIC']['fit'],
                            xic_inten=quant['XIC']['inten'],
                            as_mz=quant['apexSurvey']['mz'],
                            as_rt=quant['apexSurvey']['rt'],
                            as_least_squares=quant['apexSurvey']['fit'],
                            as_inten=quant['apexSurvey']['inten'],
                            ps_mz=quant['priorSurvey']['mz'],
                            ps_rt=quant['priorSurvey']['rt'],
                            ps_least_squares=quant['priorSurvey']['fit'],
                            ps_inten=quant['priorSurvey']['inten'],
                            pcm=quant['pcmsf'],
                            thresh=quant['thresh'],
                            score=quant['score'])

        ok = self.hdf5.appendRows(tableLoc, [quanExtraRow])

        # create combined quant data table if needed
        if self.useSecondary:
            tabName = 'quancomp'
            tableDef = 'QuanComparison'
            tableLoc = '/%s/%s' % (tabGroup, tabName)

            if not self.hdf5.hasTable(tableLoc):
                ok = self.hdf5.createTable(tabGroup, tabName, tableDef)

            quanCompRow = dict(spec_id=quant['spec_id'],
                               isolabel_id=quant['isolabel_id'],
                               pcm=quant['pcmsf'],
                               score=quant['score'],
                               primary_inten=quant['inten'],
                               primary_fit_c12=quant['exact_fitc12'],
                               primary_ls=quant['fit'],
                               primary_sum_ls=quant['cl_sum_fit'],
                               secondary_inten=quant['secondary_inten'],
                               secondary_fit_c12=quant['secondary_fitc12'],
                               secondary_ls=quant['secondary_fit'],
                               secondary_sum_ls=quant['cl_sum_fit_secondary'])

            ok = self.hdf5.appendRows(tableLoc, [quanCompRow])


    def writeQuantError(self, missingIons, IDdata):
        """
        writes missing quant ion data to the quantError table
        @param missingIons:
        @param IDdata:
        @return:
        """

        tabGroup = self.importGroup
        # tabGroup = 'rawdata'
        tabName = 'quanerror'
        tableLoc = '/%s/%s' % (tabGroup, tabName)
        tabDef = 'QuanError'

        if not self.hdf5.hasTable(tableLoc):
            ok = self.hdf5.createTable(tabGroup, tabName, tabDef)

        isos = {0: 'C12', 1: 'C13'}

        report = set(missingIons['apex'].keys())
        report.update(missingIons['prior'].keys())

        rowsToAdd = []
        colsUsed = {}

        for iso in report:
            row = dict()
            row['spec_id'] = IDdata['spec_id']
            colsUsed['spec_id'] = 1
            row['isolabel_id'] = IDdata['isolabel_id']
            colsUsed['isolabel_id'] = 1
            row['source'] = isos[iso]
            colsUsed['source'] = 1
            row['specID_apex'] = IDdata['apex']
            colsUsed['specID_apex'] = 1
            row['specID_prior'] = IDdata['prior']
            colsUsed['specID_prior'] = 1
            if iso in missingIons['apex']:
                row['foundMZ_apex'] = missingIons['apex'][iso][0]
                colsUsed['foundMZ_apex'] = 1
                row['foundInten_apex'] = missingIons['apex'][iso][1]
                colsUsed['foundInten_apex'] = 1
                row['foundMZDelta_apex'] = missingIons['apex'][iso][2]
                colsUsed['foundMZDelta_apex'] = 1
            if iso in missingIons['prior']:
                row['foundMZ_prior'] = missingIons['prior'][iso][0]
                colsUsed['foundMZ_prior'] = 1
                row['foundInten_prior'] = missingIons['prior'][iso][1]
                colsUsed['foundInten_prior'] = 1
                row['foundMZDelta_prior'] = missingIons['prior'][iso][2]
                colsUsed['foundMZDelta_prior'] = 1

            rowsToAdd.append(row)

        ok = self.hdf5.appendRows(tableLoc, rowsToAdd, colsUsed)

        pass

    def showSummary(self, allPCMs, pcmLenghts, numLabels, IsoCollect, numPeptides, allSeqs, qcSeqs, failures,
                    noMS1quant):
        """
        prints out summary of the analysis to screen
        @param allPCMs <dictionary>: containing all pcm data
        @param pcmLenghts:
        @param numLabels <dictionary>: containing the number of times each label has been linked to an XIC cluster
        @param IsoCollect <IsotopeCollection object>: containing all the isotope data for all pcms
        @param numPeptides <integer>: the number of rank 1 peptides with a score above threshold
        @param allSeqs <set>: containing the unique sequences
        @param qcSeqs <set>: containing the unique sequences with modifications to allow MS1 quantification
        @param failures <integer>: Number of peptides that have MS1 quant mods but fail quantification QC
        @param noMS1quant <integer>: Number of peptides that have no MS1 quantification mods
        @return:
        """
        print 'leastIsotopes_filter %.1f/%.1f\n' % (self.cfg.parameters['deisotoping']['maxxicls'],
                                                    self.cfg.parameters['deisotoping']['maxspecls'])
#        allPCMs = self.allPCMs
#        pcmLenghts = self.pcmLenghts
#        numLabels = self.numLabels

        freq = {}
        for pcm in allPCMs:
            num = len(allPCMs[pcm])
            if num in freq:
                freq[num] += 1
            else:
                freq[num] = 1
            isotopes = IsoCollect.getMS1quantModifiedIsoPatterns(allPCMs[pcm][0])

        for pr in [('peptides', numPeptides), ('unique seqs', len(allSeqs)), ('PCMs', len(allPCMs)),
                   ('Unique PCMseq', len(qcSeqs)), ('failures', failures), ('no MS1Quant', noMS1quant)]:
            print '#%-14s= %i' % (pr[0], pr[1])
        print '\nState\tFrequency'

        for s in isotopes['byIncreasingMZ']:
            print '%s\t%i' % (s, numLabels[s])
        print '\n#PCM\tfrequency'
        for num in freq:
            print '%3i  \t%6i' % (num, freq[num])
            #        print '\nlength\tfrequency'
            #        for l in pcmLenghts:
            #            print '%3i  \t%6i' % (l, pcmLenghts[l])

    def PCMoutput(self, pcm, quant, bestMethod, isotopes, switched, missingIons, specIDs, label2id,
                  fsFitting, fsSeparate, fsTogether):
        """
        @brief controls the output of data for each PCM
        @param pcm <PCM object>: containing the PCM quantification data
        @param quant <dictionary>: containing the quantification data
        @param bestMethod <string>: giving the method to use for quantification (XIC, apexSurvey, priorSurvey)
        @param isotopes <list>: the order the isotopes should be processed
        @param separate <optional file object>: filepath for the output file with each label on separate lines
        @param together <optional file object>: filepath for the output file with all label on teh same line
        @return:
        """

#        hdf = self.hdf
#        label2id = self.label2id
        methods = ['XIC', 'apexSurvey', 'priorSurvey']

        fileBaseData = '%s\t%s_%i\t%s\t%i\t%f' % (pcm.PCMsf, pcm.sequence, pcm.charge, bestMethod, pcm.spec_id,
                                                  pcm.score)

        fitBaseData = '%s\t%s\t%i\t|%s|' % (fileBaseData, pcm.sequence, pcm.charge, pcm.nonMS1QuantMods)

        if not quant[bestMethod]:
            return 0

        pcmData = fileBaseData
        for label in isotopes['byIncreasingMZ']:
            # step through each label
            hdf5data = dict(spec_id=pcm.spec_id, survey_id=pcm.survey, pcm=pcm.PCM, pcmsf=pcm.PCMsf, thresh=pcm.thresh,
                            score=pcm.score, isolabel_id=0, area=0.0, inten=0.0, mzdiff=0.0, ppm=0.0, monoMZ=0.0,
                            foundRT=0.0, fit=0.0, sum_theor_int=0.0, cl_sum_fit=0.0, cl_max_fit=0.0)

            if self.useSecondary:
                hdf5data['secondary_fit'] = -0.2
            else:
                hdf5data['secondary_fit'] = None

            sumIntensity = isotopes[label]['primary_sumTheo']

            # if the bestMethod has data then add to the hdf5 file
            if quant[bestMethod][label]:
                q = quant[bestMethod][label]
                hdf5data['isolabel_id'] = label2id[label]
                hdf5data['inten'] = sumIntensity * q['fitInten']
                hdf5data['corrected'] = hdf5data['inten']
                hdf5data['fit'] = q['leastSquares']
                if bestMethod == 'XIC':
                    hdf5data['mzdiff'] = q['cluster'][0]['peakmz'] - isotopes[label]['monoMZ']
                    hdf5data['foundRT'] = q['cluster'][0]['rt']
                else:
                    hdf5data['mzdiff'] = q['cluster']['specMZ'] - isotopes[label]['monoMZ']
                    hdf5data['foundRT'] = q['cluster']['specRT']

                hdf5data['ppm'] = hdf5data['mzdiff'] / isotopes[label]['monoMZ'] * 1e6

                hdf5data['cl_sum_fit'] = quant[bestMethod]['patternSetFitQuality']
                hdf5data['cl_max_fit'] = quant[bestMethod]['clusterFitMax']
                hdf5data['monoMZ'] = isotopes[label]['monoMZ']
                hdf5data['sum_theor_int'] = sumIntensity
                # hdf5data['averagine_fit'] = q['averagine']

                if self.useSecondary:

                    hdf5data['cl_sum_fit_secondary'] = quant[bestMethod]['secondaryFitQuality']
                    hdf5data['cl_max_fit_secondary'] = quant[bestMethod]['secondaryFitMax']
                    hdf5data['monoMZ'] = isotopes[label]['monoMZ']
                    hdf5data['secondary_fit'] = q['secondary_fit']
                    hdf5data['secondary_fitc12'] = q['secondary_inten']
                    hdf5data['secondary_inten'] = q['secondary_inten'] * isotopes[label]['secondary_sumTheo']
                    hdf5data['exact_fitc12'] = q['fitInten']

                if '-1ratio' in quant['XIC'][label]:
                    hdf5data['-1ratio_xic'] = quant['XIC'][label]['-1ratio']
                else:
                    hdf5data['-1ratio_xic'] = -1

                if '-1ratio' in quant['apexSurvey'][label]:
                    hdf5data['-1ratio_apex'] = quant['apexSurvey'][label]['-1ratio']
                else:
                    hdf5data['-1ratio_apex'] = -1

                if '-1ratio' in quant['priorSurvey'][label]:
                    hdf5data['-1ratio_prior'] = quant['priorSurvey'][label]['-1ratio']
                else:
                    hdf5data['-1ratio_prior'] = -1

                if bestMethod is 'XIC':
                    hdf5data['-1ratio_final'] = hdf5data['-1ratio_xic']
                elif bestMethod is 'apexSurvey':
                    hdf5data['-1ratio_final'] = hdf5data['-1ratio_apex']
                elif bestMethod is 'priorSurvey':
                    hdf5data['-1ratio_final'] = hdf5data['-1ratio_prior']
                else:
                    hdf5data['-1ratio_final'] = -2

                #                if bestMethod in ["apexSurvey", "priorSurvey"]:
                #                    hdf5data['ms1SpecID'] = quant[bestMethod][label]['ms1SpecID']
                #                else:
                #                    hdf5data['ms1SpecID'] = -1

                if bestMethod == 'apexSurvey':
                    hdf5data['ms1SpecID'] = specIDs['apex']
                elif bestMethod == 'priorSurvey':
                    hdf5data['ms1SpecID'] = specIDs['prior']
                elif specIDs['apex']:
                    hdf5data['ms1SpecID'] = specIDs['apex']
                elif specIDs['prior']:
                    hdf5data['ms1SpecID'] = specIDs['prior']
                else:
                    hdf5data['ms1SpecID'] = -1

                if bestMethod == "XIC":
                    hdf5data['source'] = bestMethod
                else:
                    hdf5data['source'] = bestMethod[:-6]

                for meth in methods:
                    if 'leastSquares' in quant[meth][label]:
                        qml = quant[meth][label]
                        if meth == 'XIC':
                            monoMZ = qml['cluster'][0]['peakmz']
                            monoRT = qml['cluster'][0]['rt']
                        else:
                            monoMZ = qml['cluster']['specMZ']
                            monoRT = qml['cluster']['specRT']

                        hdf5data[meth] = dict(fit=qml['leastSquares'], inten=qml['fitInten'] * sumIntensity, mz=monoMZ,
                                              rt=monoRT)
                    else:
                        hdf5data[meth] = dict(fit=-0.2, inten=0.0, mz=0.0, rt=0.0)

                self.writeQuant(hdf5data)

            else:
                for meth in methods:
                    hdf5data[meth] = dict(fit=-0.2, inten=0.0, mz=0.0, rt=0.0)

            qError = missingIons[label]
            specIDs['isolabel_id'] = label2id[label]
            specIDs['spec_id'] = pcm.spec_id
            if qError['apex'] or qError['prior']:
                self.writeQuantError(qError, specIDs)

            # step through the cluster identification methods and accumulate the data for output
            labelData = '%s\t%f\t%f' % (label, isotopes[label]['monoMZ'], sumIntensity)

            # now add the results from all methods
            for meth in methods:
                if quant[meth][label]:
                    q = quant[meth][label]
                    if meth == 'XIC':
                        labelData += '\t%f\t%f\t%f\t%f' % (q['cluster'][0]['peakmz'], q['cluster'][0]['rt'],
                                                           q['leastSquares'], q['fitInten'])
                    else:
                        labelData += '\t%f\t%f\t%f\t%f' % (q['cluster']['specMZ'], q['cluster']['specRT'],
                                                           q['leastSquares'], q['fitInten'])
                else:
                    labelData += '\t0.0\t0.0\t0.0\t0.0'

            pcmData += '\t%s' % labelData

            if fsFitting:
                # do expanded fitting report
                if self.useSecondary:
                    secondaryOutstring = '%.4f' % hdf5data['secondary_fit']
                else:
                    secondaryOutstring = 'switched off'

                strfit = '\t%f\t%f\t%f\t%s\t%f\t%s' % (pcm.mz, hdf5data['priorSurvey']['rt'],
                                                       hdf5data['apexSurvey']['rt'], label, hdf5data['foundRT'],
                                                       secondaryOutstring)

                sumTheo = 0
                theoMZ_outStr = ''
                theoRelInts_outStr = ''

                # add the theoretical mass and intensity data
                for repIso in self.cfg.parameters['output']['reportisos']:
                    if repIso in isotopes[label]['MZs']:
                        theoMZ_outStr += '\t%f' % (isotopes[label]['MZs'][repIso])
                        theoRelInts_outStr += '\t%f' % (isotopes[label]['primary'][repIso])
                    else:
                        theoMZ_outStr += '\t'
                        theoRelInts_outStr += '\t'
                strfit2 = '\t%f\t%f\t%f\t%i' % (isotopes[label]['monoMZ'], hdf5data['monoMZ'], hdf5data['inten'],
                                                switched)

                methData = ''
                for src in methods:
                    usedInts_outStr = ''

                    if quant[src][label]:
                        leastsq = quant[src][label]['leastSquares']
                        inten = quant[src][label]['fitInten'] * sumIntensity
                    else:
                        leastsq = -0.2
                        inten = 0

                    for repIso in self.cfg.parameters['output']['reportisos']:
                        if quant[src][label] and repIso in quant[src][label]['cluster']:
                            usedInts_outStr += '\t%f' % quant[src][label]['cluster']['usedInten'][repIso]
                        else:
                            usedInts_outStr += '\t0'

                    methData += '\t%f%s\t%f' % (leastsq, usedInts_outStr, inten)

                outstring = '%s%s%s%s%s%s\n' % (fitBaseData, strfit, theoMZ_outStr, theoRelInts_outStr, strfit2,
                                                methData)
                self.fsFitting.write(outstring)
                pass

            if fsSeparate:
                # do simplified report with 1 line per PCM/label
                self.fsSeparate.write('%s\t%s\n' % (fileBaseData, labelData))
        if fsTogether:
            # do simplified report with all labels on one line for each PCM
            self.fsTogether.write('%s\n' % pcmData)
        return

    def createOutputFiles(self, suffix):

        cfg = self.cfg
        hdfname = self.hdf5.filePath.stem
        if suffix:
            suffix = '_' + suffix

        outdir = cfg.parameters['runtime']['datadir'].joinpath(cfg.parameters['output']['outdir'])
        if not outdir.exists():
            outdir.mkdir(parents=True)

        fitting = outdir.joinpath('%s_%s%s.txt' % (hdfname, cfg.parameters['output']['fitting'], suffix))
        separate = outdir.joinpath('%s_%s%s.txt' % (hdfname, cfg.parameters['output']['separate'], suffix))
        together = outdir.joinpath('%s_%s%s.txt' % (hdfname, cfg.parameters['output']['together'], suffix))

        # open the fitting file and add header data
        # self.fsFitting = fitting.open('w')
        self.fsFitting = open(str(fitting), 'w')
        strMZ = ''
        strInt = ''
        strFound = ''
        strMatched = ''
        strUsed = ''

        str_theoMZ = ''
        str_theoInt = ''
        self.fsFitting.write('pcm\tp_c\tsource\tspec_id\tscore\tSeq\tcg\tmod\tmz\tRTmsms\tRTapex\tlabel' +
                             '\tfoundRT\taveragine_LS')

        for repIso in self.cfg.parameters['output']['reportisos']:
            if repIso < 0:
                str_theoMZ += '\tmz_%i' % (abs(repIso))
                str_theoInt += '\ttherInt_%i' % (abs(repIso))
            else:
                str_theoMZ += '\tmz%i' % repIso
                str_theoInt += '\ttherInt%i' % repIso

        self.fsFitting.write(str_theoMZ + str_theoInt + '\tmonoMZ_th\tmonoMZ_me\tfinal_quantVal\tswitched')

        for src in ('XIC', 'apexSurvey', 'priorSurvey'):
            str_usedInt = ''

            for repIso in self.cfg.parameters['output']['reportisos']:
                if repIso < 0:
                    str_usedInt += '\t%s_used_%i' % (src, abs(repIso))
                else:
                    str_usedInt += '\t%s_used%i' % (src, repIso)

            self.fsFitting.write("\t%s_LS%s\t%s_quantVal" % (src, str_usedInt, src))
        self.fsFitting.write("\n")

        # open the separate file and write headers
        # self.fsSeparate = separate.open('w')
        self.fsSeparate = open(str(separate), 'w')
        self.fsSeparate.write('pcm\tp_c\tsource\tspec_id\tscore\tlabel\ttheoryMZ\tsumTheory\tXIC_mz\tXIC_rt\tXIC_LS' +
                              '\tXIC_Inten\tAS_mz\tAS_rt\tAS_LS\tAS_Inten\tPS_mz\tPS_rt\tPS_LS\tPS_Inten\n')

        # open the together file and write headers
        # self.fsTogether = together.open('w')
        self.fsTogether = open(str(together), 'w')
        self.fsTogether.write('pcm\tp_c\tsource\tspec_id\tscore\tlabel_L\ttheoryMZ_L\tsumTheory_L\tXIC_mz_L\tXIC_rt_L' +
                              '\tXIC_LS_L\tXIC_Inten_L\tAS_mz_L\tAS_rt_L\tAS_LS_L\tAS_Inten_L\tPS_mz_L\tPS_rt_L' +
                              '\tPS_LS_L\tPS_Inten_L')
        self.fsTogether.write('\tlabel_M\ttheoryMZ_M\tsumTheory_M\tXIC_mz_M\tXIC_rt_M\tXIC_LS_M\tXIC_Inten_M\tAS_mz_M' +
                              '\tAS_rt_M\tAS_LS_M\tAS_Inten_M\tPS_mz_M\tPS_rt_M\tPS_LS_M\tPS_Inten_M')
        self.fsTogether.write('\tlabel_H\ttheoryMZ_H\tsumTheory_H\tXIC_mz_H\tXIC_rt_H\tXIC_LS_H\tXIC_Inten_H\tAS_mz_H' +
                              '\tAS_rt_H\tAS_LS_H\tAS_Inten_H\tPS_mz_H\tPS_rt_H\tPS_LS_H\tPS_Inten_H\n')

        return

    def closeFiles(self):
        self.fsFitting.close()
        self.fsSeparate.close()
        self.fsTogether.close()
