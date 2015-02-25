"""
This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This file handles the processing and storage of spectral data from Xcalibur
.raw files.


isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/cellzome/isobarquant
"""

# python modules
import os
import numpy as np
import math
import traceback
import shelve
import sys

sys.path.insert(0, '..')
# cellzome CommonUtils
from CommonUtils.XICprocessor import XIC
from CommonUtils.MathTools import Statistics
from CommonUtils.tools import *
from CommonUtils.ionmanipulation import IonProcessing
import CommonUtils.progressbar as progBar
import CommonUtils.ExceptionHandler as ExHa

# pyMSsafe modules
from SpectrumProcessor import Processor, IsotopeModel

stats = Statistics()


class Datamanager():
    def __init__(self, xraw, cfg, logs, hdf, quanmeth, tempFile):
        """
        @brief initialises the datamanager object
        @param xraw <object>: the XRawFile object - interface for the Xcalibur data
        @param cfg <object>: containing the configuration data from file (and database when sorted)
        @param hdf <object>: the interface to the hdf5 output file
        @param quanmeth <string>: the name of the quantitation method
        @param tempFile <string>: file path string for the shelve data
        """
        # self.w = WMI('.')
        self.pid = os.getpid()

        self.raw = xraw
        self.cfg = cfg
        self.logs = logs
        self.hdf = hdf
        self.quanmeth = quanmeth
        self.tempFile = tempFile
        self.specproc = Processor(cfg)
        self.specs = None
        self.mssurvey = {}
        global isomodel
        isomodel = IsotopeModel(cfg, logs, 200, 0.01, quanmeth)
        self.ishcd = 0
        self.peaks3D = {}
        self.bins = 0
        self.maxspec = 0
        self.unitnum = 0
        self.currmeth = 0
        self.offsets = []
        self.lastQuant = (0, 0)
        self.server = ''
        self.job = ''
        self.job_id = ''
        self.client = ''
        self.wf_id = ''
        self.instrument = ''
        self.surveydata = None
        self.isOrbi = True
        self.units = None
        self.probunits = None
        self.acts = None
        self.quan = 0
        self.isotopeIDs = None
        self.coalescence = None
        self.isotopes = None
        self.IonProc = IonProcessing(cfg, self.logs.log)

    def addLUXdata(self, job, pid):
        self.job = job
        self.job_id = job.job_id
        self.wf_id = 0

    def processSpectra(self, quan, watch, rawbase):
        """
        @brief controls the analysis of the raw file
        @param quan <dictionary>: containing the quantiation masses if any
        @param watch <object>: stopwatch object that recods timings for the analysis
        @param rawbase <string>: the stripped filename to use to find the fragmentationmethod
        """
        # functions
        loadspec = self.loadSpectrum
        updatemsms = self.updateMSMS
        calcbgisos = self.calcBackgroundIsotopes
        fitisotopes = self.fitIsotopes
        processheaderinfo = self.processHeaderInfo
        wfid = self.wf_id

        # data references
        maxspec = self.maxspec
        cfg = self.cfg
        xRaw = self.raw
        inst = xRaw.getInstMethodNames()[0]
        self.instrument = inst

        num_specs = self.loadParameters(quan)
        width = self.getIsolationWidth()

        # move s2idata to a shelve object
        s2idata = {}
        s2iShelf = shelve.open(self.tempFile, protocol=2)

        # needs to be referenced after setting the specs array
        msmsspecs = self.specs
        ishcd = self.ishcd

        # test the frag method in the file fits to a defined frag method
        acts = self.acts[1]
        num = len(acts)

        # get the fragmeths from the parameter file
        setacts = cfg.formatFragMeths()
        for meth in setacts:
            if len(meth) == num:
                ok = 0
                for i in range(num):
                    if meth[i][0] == acts[i]:
                        ok += 1
                if ok == num:
                    # method matches
                    break

        # test the fragmeth against the file methods
        fail = 1
        if len(meth) == num:
            ok = 0
            for i in range(num):
                if meth[i][0] == acts[i]:
                    ok += 1
            if ok == num:
                # method matches so update the units data with the spectrum usage
                for m in range(len(meth)):
                    self.units[1][m]['use'] = meth[m][1]
                fail = 0

        if fail:
            # raise an error to report to LUX
            raise ExHa.FragmentMethodError('Missmatched: Actual %s  vs Expected %s' % (str(acts), str(meth)))

        self.writeUnitData(self.units, self.probunits)

        watch.rec('Load parameters')
        if maxspec:
            self.logs.log.info('processing %d spectra' % maxspec)
        else:
            self.logs.log.info('processing %d spectra' % num_specs)

        # load the spectrum data
        specs = xRaw.genSpecNum()
        if maxspec:
            pBar = progBar.ProgressBar(widgets=progBar.name_widgets, maxval=maxspec, name='Spectra anal.').start()
        else:
            pBar = progBar.ProgressBar(widgets=progBar.name_widgets, maxval=num_specs, name='Spectra anal.').start()
        # cnt.sec()
        surveyid = 0
        s2ilate = {}
        spRow = -1
        lastUpdated = -1

        nspec = 0
        specNum = 0
        try:
            while 1:
                specNum = specs.next()
                nspec += 1
                # rtn = cnt.next()
                pBar.update(nspec)

                if nspec % 1000 == 0:
                    s2iShelf.sync()

                    # check if this is a LUX job - report to LUX when dot is written
                    if wfid:
                        self.updateLux()

                # load the spectrum and noise data
                spdata, spRow = loadspec(specNum, surveyid, spRow)
                noise = xRaw.getNoiseData(specNum)
                if len(noise) > 0:
                    self.writeNoise(specNum, noise)

                if spdata == -1:
                    # deal with skipped spectra: rejected scan event or missmatched headder data
                    continue
                elif spdata == 0:
                    # this is MS/MS data
                    pBar.nextSecondary()
                    # only calculate the s2i data once for each unit
                    unit = msmsspecs[spRow]['unit']

                    if unit in s2idata:
                        # unit located: copy data to new spectrum
                        s2idata[unit]['msms']['specs'][specNum] = dict(specs=specNum, rt=msmsspecs[spRow]['rt'],
                                                                       row=spRow, use=self.currmeth['use'])
                    else:
                        # find the data ready for s2i calculations
                        if msmsspecs[spRow]['preccharge'] == 0:
                            background = 0.0
                            foundisos = {}
                            isoions = []
                            c12 = dict(fit=0, base=0, c12=0, inten=0, mz=0, norm=0, offset=0, smooth=0)
                        else:
                            background, foundisos, isoions = calcbgisos(spRow, surveydata)
                            c12 = fitisotopes(foundisos, msmsspecs[spRow]['setmass'], msmsspecs[spRow]['preccharge'])
                        thresh = surveydata.getThreshold(msmsspecs[spRow]['setmass'])
                        s2idata[unit] = dict(setmass=msmsspecs[spRow]['setmass'], monomz=msmsspecs[spRow]['monomz'],
                                             charge=msmsspecs[spRow]['preccharge'])
                        s2idata[unit]['early'] = dict(bkgd=background, foundisos=foundisos, rt=surveydata.header['rt'],
                                                      thresh=thresh, c12=c12, surveyid=surveyid)
                        msms = {specNum: dict(specs=specNum, row=spRow, rt=msmsspecs[spRow]['rt'],
                                              use=self.currmeth['use'])}
                        s2idata[unit]['msms'] = dict(specs=msms, isoions=isoions)
                        msmsspecs[spRow]['threshearly'] = thresh

                    s2ilate[unit] = s2idata[unit]['msms']['specs'].keys()
                else:
                    # for HCD data make sure the IDspec and QuanSpec are set properly
                    spdata.noise = noise
                    if ishcd and specNum - surveyid > 1:
                        updatemsms(lastUpdated + 1, spRow + 1)
                        lastUpdated = spRow
                    surveyid = specNum
                    self.surveydata = spdata
                    surveydata = spdata

                    # append the late data for all the MS/MS spectra since the previous survey
                    for unit in s2ilate:
                        s2iunit = s2idata[unit]
                        s2ispecs = s2iunit['msms']['specs']

                        # update the quan and id spec data
                        idSpec = 0
                        quantSpec = 0
                        for sn in s2ilate[unit]:
                            specdata = s2ispecs[sn]
                            use = specdata['use']
                            sp = (specdata['row'], specdata['specs'])

                            if use == 'IQ':
                                msmsspecs[sp[0]]['quanspec'] = sp[1]
                                msmsspecs[sp[0]]['idspec'] = sp[1]
                            elif use == 'I':
                                idSpec = sp
                                msmsspecs[sp[0]]['idspec'] = sp[1]
                                if quantSpec:
                                    msmsspecs[sp[0]]['quanspec'] = quantSpec[1]
                                    msmsspecs[quantSpec[0]]['idspec'] = idSpec[1]
                            elif use == 'Q':
                                quantSpec = sp
                                msmsspecs[sp[0]]['quanspec'] = sp[1]
                                if idSpec:
                                    msmsspecs[sp[0]]['idspec'] = idSpec[1]
                                    msmsspecs[idSpec[0]]['quanspec'] = quantSpec[1]

                        # find the data ready for s2i calculations
                        thresh = surveydata.getThreshold(s2iunit['setmass'])
                        if unit == 13:
                            pass
                        if s2iunit['charge'] == 0:
                            c12 = dict(fit=0, base=0, c12=0, inten=0, mz=0, norm=0, offset=0, smooth=0)
                            background = 0.0
                            foundisos = {}
                        else:
                            row = s2ispecs[s2ilate[unit][0]]['row']
                            background, foundisos, isoions = calcbgisos(row, surveydata)
                            c12 = fitisotopes(foundisos, msmsspecs[row]['setmass'], msmsspecs[row]['preccharge'])

                        s2iunit['late'] = dict(bkgd=background, foundisos=foundisos, rt=surveydata.header['rt'],
                                               thresh=thresh, c12=c12, surveyid=surveyid)
                        processheaderinfo(s2iunit)

                        # now shelve the unit and take it out of the s2idata dictionary
                        s2iShelf[str(unit)] = s2idata.pop(unit)

                    s2ilate = {}

                if maxspec and specNum > maxspec:
                    raise StopIteration
        except StopIteration:
            if s2ilate:
                # deal with the possibility that there was no last survey
                for unit in s2ilate:
                    s2iunit = s2idata[unit]
                    s2iunit['late'] = s2iunit['early'].copy()
                    s2iunit['late']['rt'] = s2iunit['msms']['specs'][max(s2iunit['msms']['specs'])]['rt'] + 1.0
                    processheaderinfo(s2iunit)
            # cnt.end()
            pBar.finish()
            self.logs.log.info('Creating HDF5 indexes')
            self.createIndexes()
            if wfid:
                self.updateLux()

            self.logs.log.info('spectrum processing finished')
            xRaw.close()
        except:
            xRaw.close()
            trace = repr(traceback.format_exception(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
            return {'code': 2, 'error': 'pyMSsafe Spectrum processing error (spec %d): %s' % (specNum, trace)}
        if self.offsets:
            msd = stats.meanStd(self.offsets)
            if msd[1] is None:
                self.logs.log.info('average Reporter Ion offset = %.2f mDa +/- None' % (msd[0] * 1000))
            else:
                self.logs.log.info('average Reporter Ion offset = %.2f mDa +/- %.2f' % (msd[0] * 1000, msd[1] * 1000))
        watch.rec('Process Spectra')

        self.specs = msmsspecs
        s2iShelf.close()

        return {'code': 0}

    def getIsolationWidth(self):
        """
        @brief gets the isolation width data from the Unit data and interpolates real values from the config data
        @return:
        """
        units = self.units
        s2iParam = self.cfg.parameters['sig2interf']
        isolation = s2iParam['isolation']

        widths = {}
        for u in range(1, len(units)):
            for frag in units[u]:
                widths.setdefault(frag['activation'], []).append(frag['isolation'])

        for k in widths:
            widths[k] = min(widths[k])

        # use the hcd width if present otherwise use the narrowest value
        if 'HCD' in widths:
            width = widths['HCD']
        else:
            width = min(widths.values())

        if width in isolation:
            values = isolation[width]
            s2iParam['winfull'] = values[0]
            s2iParam['winscaled'] = values[1]
        else:
            ordered = isolation.keys()
            ordered.sort()

            msg = 'Isolation width %.1f beyond %s range, extrapolating from %.1f'
            if width < ordered[0]:
                # width too low: extrapolate
                idx = 0
                isol1 = ordered[idx]
                self.logs.log.info(msg % (width, 'min', ordered[0]))
            elif width > ordered[-1]:
                # width too high: extrapolate
                idx = len(ordered) - 2
                isol1 = ordered[idx]
                self.logs.log.info(msg % (width, 'max', ordered[-1]))
            else:
                for idx, isol1 in enumerate(ordered):
                    if isol1 < width < ordered[idx + 1]:
                        break

            scaled1 = isolation[isol1][0]
            full1 = isolation[isol1][1]
            isol2 = ordered[idx + 1]
            scaled2 = isolation[isol2][0]
            full2 = isolation[isol2][1]

            s2iParam['winfull'] = (full2 - full1) / (isol2 - isol1) * (width - isol1) + full1
            s2iParam['winscaled'] = (scaled2 - scaled1) / (isol2 - isol1) * (width - isol1) + scaled1

        return width

    def processHeaderInfo(self, headerinfo):
        """
        """
        # functions
        interpolate = self.interpolate

        # data
        msmsspecs = self.specs
        headlate = headerinfo['late']['c12']
        headearly = headerinfo['early']['c12']
        mzearly = headearly['mz']
        mzlate = headlate['mz']
        specs = headerinfo['msms']['specs']
        params = self.cfg.parameters
        neutron = float(params['general']['neutron'])

        xicmin = 10000
        xicmax = 0
        sumrt = 0.0
        for row in specs:
            specdata = specs[row]
            sumrt += specdata['rt']
            if specdata['rt'] < xicmin:
                xicmin = specdata['rt']
            if specdata['rt'] > xicmax:
                xicmax = specdata['rt']

            specs[row]['thresh'] = interpolate(headerinfo['early']['thresh'], headerinfo['early']['rt'],
                                               headerinfo['late']['thresh'], headerinfo['late']['rt'],
                                               specdata['rt'])
        headerinfo['rt'] = sumrt / len(specs)
        headerinfo['xicmin'] = xicmin - float(params['xic']['afterpeak'])
        headerinfo['xicmax'] = xicmax + float(params['xic']['beforepeak'])

        delta = abs(mzlate - mzearly)
        if delta > 0.05:
            # 12C isos have different masses
            if mzearly == 0:
                # no early data so use late
                precmz = mzlate
                offset = headlate['offset']
                norm = headlate['norm']
                maxintc13 = headlate['base']
                surv_fit = headlate['fit']
            elif mzlate == 0:
                # no late data so use early
                precmz = mzearly
                offset = headearly['offset']
                norm = headearly['norm']
                maxintc13 = headearly['base']
                surv_fit = headearly['fit']
            elif headearly['fit'] <= headlate['fit']:
                # early data has better isotopic fit
                precmz = mzearly
                offset = headearly['offset']
                norm = headearly['norm']
                maxintc13 = headearly['base']
                surv_fit = headearly['fit']
            elif headearly['fit'] > headlate['fit']:
                # late data has better isotopic fit
                precmz = mzlate
                offset = headlate['offset']
                norm = headlate['norm']
                maxintc13 = headlate['base']
                surv_fit = headlate['fit']
        elif mzearly == 0:
            # masses are both zero
            precmz = 0
            offset = 0
            norm = 0
            maxintc13 = headerinfo['monomz']
            surv_fit = 0
        else:
            # masses are the same use the best fitting
            if headearly['fit'] <= headlate['fit']:
                # early data has better isotopic fit
                precmz = mzearly
                offset = headearly['offset']
                norm = headearly['norm']
                maxintc13 = headearly['base']
                surv_fit = headearly['fit']
            elif headearly['fit'] > headlate['fit']:
                # late data has better isotopic fit
                precmz = mzlate
                offset = headlate['offset']
                norm = headlate['norm']
                maxintc13 = headlate['base']
                surv_fit = headlate['fit']

        for sp in specs:
            msmsspecs[specs[sp]['row']]['precmz_surv'] = precmz
        headerinfo['offset_surv'] = offset
        headerinfo['norm_surv'] = norm
        headerinfo['precmz_surv'] = precmz
        headerinfo['surv_fit'] = surv_fit

        # adjust the isotopes to the precmz_surv
        if headerinfo['precmz_surv'] > 0:
            for i in range(len(headerinfo['msms']['isoions'])):
                iso = headerinfo['msms']['isoions'][i][0]
                headerinfo['msms']['isoions'][i] = (iso, precmz + (iso - offset) * neutron / headerinfo['charge'])

            for iso in headerinfo['msms']['isoions']:
                if iso[0] == maxintc13:
                    headerinfo['basemass'] = iso[1]
                    headerinfo['baseoffset'] = iso[0]
        else:
            headerinfo['basemass'] = headerinfo['monomz']
            headerinfo['baseoffset'] = 0

    @staticmethod
    def interpolate(valueearly, timeearly, valuelate, timelate, timeinterp):
        """
        @brief interpolates a value at a time from two values at time points either side of the desired time
        @param valueearly <float>: value at the early time point
        @param timeearly <float>: time for the early time point
        @param valuelate <float>: value at the late time point
        @param timelate <float>: time for the late time point
        @param timeinterp <float>: time for the interpolation point
        @return interp <float>: the interpolated value
        """

        interp = (timeinterp - timeearly) / (timelate - timeearly) * (valuelate - valueearly) + valueearly
        return interp

    def fitIsotopes(self, isotopes, mz, charge):
        """
        @brief takes isotope data optimises the cluster fitting
        @param isotopes <list>: containing all the isotope data (offset, intensity etc)
        @param mz <float>: mass to use to get cluster intensity ratios
        @param charge <float>: charge state of the cluster
        @return c12 <tuple>: containing the optimum 12C ion from the isotope data
        """
        nomatch = dict(c12=0.0, smooth=0.0, inten=0.0, offset=0,
                       norm=0.0, fit=0.0, mz=0.0)
        if not isotopes:
            # no data to work with
            return nomatch

        isoset = isomodel.findIsotopeSet(mz * charge, 2)
        offsets = isotopes.keys()
        offsets.sort()
        deisoMeth = self.cfg.parameters['deisotoping']['deisotope']

        # check the the possibility of offset isotope pattern matching
        if deisoMeth == 'leastsq':
            fit = self.isoFitLeastSQ(isoset['inten'], isotopes, charge, 1)
            if fit:
                if fit[3] in isotopes:
                    best = isotopes[fit[3]]
                    ion = dict(c12=best[0], smooth=best[1], inten=best[2], offset=best[3],
                               norm=fit[2], fit=fit[0], mz=fit[4], base=isoset['max'])
                else:
                    return nomatch
            else:
                # no best fit data
                return nomatch
        elif deisoMeth == 'absdelta':
            best = self.isoFitAbsoluteDelta(isoset['inten'], isotopes)
            ion = dict(c12=best[0], smooth=best[1], inten=best[2], offset=best[3], norm=best[2], fit=0,
                       base=isoset['max'])
        else:
            raise ExHa.MSdataConsistancyError('No isotope fitting method selected in config file.')

        return ion

    def isoFitAbsoluteDelta(self, isos, ions):
        logtol = math.log10(self.cfg.parameters['deisotoping']['low_int_tol'])
        best = []
        isokeys = ions.keys()
        isokeys.sort()
        ioff = isokeys.pop(0)

        while ioff < 1:
            # test combinations starting from the lowest mass
            c12 = ions[ioff]
            run = 0
            avratio = 0.0
            for i in range(min(len(isokeys), len(isos))):
                ratio = math.log10(ions[isokeys[i]][1] / c12[1] / isos[i])
                if abs(ratio) <= logtol:
                    run += 1
                    avratio += ratio
                else:
                    break
            if run > 0:
                avratio = abs(avratio) / run

            if not best:
                best = [(run, avratio, ioff)]
            elif run > best[0][0]:
                best = [(run, avratio, ioff)]
            elif run == best[0][0]:
                best.append((run, avratio, ioff))

            if isokeys:
                ioff = isokeys.pop(0)
            else:
                break

        if best:
            best.sort()
            if len(best) > 1 and best[0][0] == 0:
                # if no 13C data then choose the most intense ion
                c12 = ions[best[0][2]]
                for b in best:
                    if ions[b[2]][1] > c12[1]:
                        c12 = ions[b[2]]
            else:
                # otherwise choose the best fitting (lowest avratio)
                b = best[0]
                c12 = ions[b[2]]
        else:
            c12 = (0.0, 0.0, 0.0)

        return c12

    def isoFitLeastSQ(self, isos, ions, charge, havestart):
        """
        @brief find the best isotope fit to the ion data using least squares optimisation
        @param isos <list>: of the relative intensities of the theoretical averagine isotope model
        @param ions <dictionary: containing the observed ion data indexed by the isotope offset
        @param charge <integer>: charge state of the ion to be deisotoped
        @param havestart <integer>: flag indicating if the data should have the start ion present
        @return ion <tuple>: containing the mz, intensity and offset of the 12C ion of the prefered fit
        """
        calcSumSQ = self.calcSumSQ
        selectLowest = self.selectLowest
        neutron = self.cfg.parameters['general']['neutron']
        lsLimit = self.cfg.parameters['deisotoping']['ls_limit']
        intLimit = self.cfg.parameters['deisotoping']['int_limit']

        num = len(isos)
        fits = []

        # start points limited to clusters that could contain the monoisotopic ion
        for start in range(0, 1 - num, -1):
            # test all start points if they have ions corresponding to the start
            if havestart and start not in ions:
                break
            intens = []
            minint = 1e10
            maxint = 0
            matches = []
            summass = 0.0
            summassint = 0.0
            sumint = 0.0
            for i in range(num):
                if start + i in ions:
                    newint = ions[start + i][1]
                    matches.append(1)
                    summass += (ions[start + i][0] - i * neutron / charge)
                    summassint += (ions[start + i][0] - i * neutron / charge) * ions[start + i][1]
                    sumint += ions[start + i][1]
                else:
                    newint = 1
                    matches.append(0)
                if newint > maxint:
                    maxint = newint
                if newint < minint:
                    minint = newint
                intens.append(newint)
            if maxint > 0:
                opts = []
                nextPos = -1
                opts.append((minint, calcSumSQ(isos, intens, minint)))
                midint = (maxint + minint) / 2
                opts.append((midint, calcSumSQ(isos, intens, midint)))
                opts.append((maxint, calcSumSQ(isos, intens, maxint)))

                while opts[0][1] < opts[1][1]:
                    # first point is lowest - need another point lower
                    newint = opts[0][0] / 2
                    opts.append((newint, calcSumSQ(isos, intens, newint)))
                    opts.sort()
                    opts, nextPos = selectLowest(opts)

                while opts[2][1] < opts[1][1]:
                    # last point is lowest - need another point
                    newint = opts[2][0] * 2
                    opts.append((newint, calcSumSQ(isos, intens, newint)))
                    opts.sort()
                    opts, nextPos = selectLowest(opts)

                if nextPos == -1:
                    opts, nextPos = selectLowest(opts)

                while (opts[nextPos][1] - opts[1][1]) / opts[1][1] > lsLimit and \
                      (opts[2][0] - opts[0][0]) / opts[1][0] > intLimit:
                    midint = (opts[1][0] + opts[nextPos][0]) / 2
                    opts.append((midint, self.calcSumSQ(isos, intens, midint)))
                    opts.sort()
                    opts, nextPos = selectLowest(opts)
                fits.append((opts[1][1], matches, opts[1][0], start, summassint / sumint))
        if fits:
            fits.sort()
            return fits[0]
        else:
            return 0, [0, 0], 0, 0, 0

    @staticmethod
    def selectLowest(data):
        if len(data) < 4:
            pass
        elif data[1][1] < data[2][1]:
            data = data[:3]
        else:
            data = data[1:]

        if data[0][1] > data[2][1]:
            nextPos = 0
        else:
            nextPos = 2
        return data, nextPos

    @staticmethod
    def calcSumSQ(isos, ions, norm):

        sumsq = 0.0
        for i in range(len(ions)):
            sumsq += (isos[i] - ions[i] / norm) ** 2

        return sumsq

    def loadParameters(self, quantMethod):
        """
        @brief loads the parameter information from the raw file and transfers to the hdf5 file
        @param quan <dictionary>: containing the quantiation masses if any
        @return length <integer>, first <integer>, last <integer>: the number of spectra and the number of the first
        and last spectrum
        """
        cfg = self.cfg

        wanted = cfg.parameters['wanted']

        xRaw = self.raw

        if self.instrument == 'LTQ':
            self.isOrbi = True
        else:
            self.isOrbi = False

        # HPLC method data
        hplc = xRaw.getLCMethods(cfg.parameters['methodnames'])
        if hplc['analytical']:
            self.writeLCparameters(hplc, wanted)

        # MS method data
        ms, order = self.raw.getMSMethods(cfg.parameters['methodnames']['ms'])

        ishcd = []
        for ev in ms['activation']:
            for act in ms['activation'][ev]:
                if 'HCD' in act:
                    ishcd.append(len(ms['activation'][ev]))
                else:
                    ishcd.append(0)
        self.ishcd = max(ishcd)
        self.units = ms.pop('units')
        if 'units' in order:
            order.remove('units')
        self.acts = ms['activation'].copy()
        self.probunits = ms.pop('unit problems')
        if 'units problems' in order:
            order.remove('unit problems')

        # build list of scan events
        scanevents = []
        for i in range(ms['num scan events']):
            scanevents.append('Scan Event %d' % (i + 1))

        if self.isOrbi:
            subsets, scanevents = self.fixScanEvents(wanted['ms_subset'], scanevents)
            self.writeMSparameters(ms, wanted['ms_param'], subsets, order)
        else:
            subsets, scanevents = self.fixScanEvents(wanted['exactive_subset'], scanevents)
            self.writeExactiveMSparameters(ms, wanted['exactive_param'], subsets, order)

        self.logs.log.info('MS has %d parameters' % (len(ms)))

        # Tune data
        try:
            tune = xRaw.getTuneParam()
            self.writeTuneParameters(tune, cfg.parameters['wanted']['tune_param'],
                                     cfg.parameters['wanted']['tune_subset'])
            self.logs.log.info('Tune has %d parameters' % (len(tune)))
        except:
            self.logs.log.info('No tuning parameters for file')
        # else: print 'No tuning parameters for file'

        numSpecs, first, last = xRaw.getNumSpectra()

        spGen = xRaw.genSpecNum()
        length = 0
        try:
            while 1:
                sp = spGen.next()
                filt = xRaw.getFilterForScanNum(sp)
                if filt.find(' ms ') == -1:
                    length += 1
        except StopIteration:
            # finish
            pass

        columns = [('spec', int), ('unit', int), ('order', int), ('scanevent', int), ('monomz', float),
                   ('precmz', float), ('precmz_surv', float), ('preccharge', int),
                   ('precsd', float), ('setmass', float), ('survey', int), ('quanspec', int),
                   ('idspec', int), ('rt', float), ('scanapex', int), ('rtapex', float),
                   ('inten', float), ('area', int), ('fwhm', float), ('c13iso', int),
                   # ('s2iearly', float), ('s2ilate', float), ('s2i', float),
                   # ('c12early', float), ('c12late', float), ('c12', float),
                   ('mzearly', float), ('mzlate', float),
                   ('threshearly', float), ('threshlate', float), ('thresh', float),
                   ('fragmeth', 'S5'), ('fragenergy', float), ('sumrepint', float), ('sumreparea', float)]

        self.specs = np.ndarray(length, dtype=columns)
        for r in range(length):
            self.specs[r]['spec'] = 0

        # sort out the quan masses if any
        self.quan = 0
        if quantMethod:
            if quantMethod['source'] == 'ms1':
                self.quan = 0
            else:
                self.quan = quantMethod

            # have all the quant parameters locally
            quantTol = cfg.parameters['quantitation']

            if self.ishcd:
                # HCD data so set tolerences for HCD accuracy
                quantMethod.update(dict(search=quantTol['hcd_wide_tol'], match=quantTol['hcd_tol']))
            else:
                # non-HCD low accuracy data
                quantMethod.update(dict(search=quantTol['trap_wide_tol'], match=quantTol['trap_tol']))
            self.writeIsotopes(quantMethod)

            # extract isotope IDs and any potential coalescnce pairs
            idList = []
            isotopes = quantMethod['quantmasses']
            for key in isotopes:
                idList.append(key)

            idList.sort(key=lambda x: isotopes[x][0]['mass'])
            coalescence = []
            i = 0
            while i < len(idList) - 1:
                thisID = idList[i]
                nextID = idList[i + 1]

                if abs(isotopes[thisID][0]['mass'] - isotopes[nextID][0]['mass']) < 0.10:
                    # ions within 10 mDa could be subject to coalescence
                    midMass = (isotopes[thisID][0]['mass'] + isotopes[nextID][0]['mass']) / 2
                    coalescence.append((thisID, nextID, midMass))
                    i += 1

                i += 1

            self.isotopeIDs = idList
            self.coalescence = coalescence

            if quantMethod['source'] == 'ms2':
                # this is MS2 quant and only one mass per isotope
                for key in idList:
                    if len(isotopes[key]) == 1:
                        isotopes[key] = isotopes[key][0]
                    else:
                        raise ExHa.MSdataConsistancyError('too many items in quant label')

            self.isotopes = quantMethod.copy()

        return numSpecs

    @staticmethod
    def fixScanEvents(subsets, scanevents):
        if 'Scan Events' in subsets:
            subsets.remove('Scan Events')
            subsets += scanevents
        return subsets, scanevents

    def writeLCparameters(self, hplc, wanted):
        """
        @brief organises the LC parameter data for writing to HDF5 file
        @param hplc <dictionary>: containing the parameters from the LC
        @param wanted <dictionary>: containing all the instrument parameters to be included
        @return:
        """

        writing = []
        profile = []
        if hplc['lctype'] == 'Eksigent':
            # LTQOrbi uses Eksigent Nano LC
            profileTable = 'Profile'
            wantLC = wanted['lc_param']
            for section in ['lctype', 'analytical', 'loading']:
                if section == 'lctype':
                    writing.append(dict(set='Eksigent', subset='main', parameter=section, value=hplc[section]))
                else:
                    mainset = 'LC %s' % section
                    HPLCdata = hplc[section]
                    for key in HPLCdata:
                        if key in ['Qa_Profile', 'Qb_Profile']:
                            # these are LC conditions
                            writing.append(dict(set=mainset, subset=key, parameter='numpts',
                                                value=HPLCdata[key]['numpts']))
                            profid = mainset[3:7] + key[:2]
                            writing.append(dict(set=mainset, subset=key, parameter='profile_id', value=profid))

                            for gdt in HPLCdata[key]['pts']:
                                profile.append(dict(id=profid, time=gdt['t'], time_unit=gdt['t_unit'], flow=gdt['flow'],
                                                    flow_unit=gdt['flow_unit']))
                        elif key in wantLC:
                            writing.append(dict(set=mainset, subset='main', parameter=key, value=HPLCdata[key]))

                    self.logs.log.info('LC-%s has %d parameters' % (section, len(HPLCdata)))

        elif hplc['lctype'] == 'Dionex':
            # Exactive uses Dionex LC
            profileTable = 'Dionex'
            wantDionex = wanted['dionex_param']

            for key in hplc:
                if key in ['analytical', 'loading']:
                    # these are the pump specific data
                    for param in hplc[key]:
                        if param == 'profile':
                            # this is the profile data
                            lcID = 'DionexLC' + key[:4]
                            for point in hplc[key]['profile']:
                                profile.append(dict(id=lcID, time=point['time'], time_unit='min', percent_b=point['b'],
                                                    flow=point['flow'], flow_unit=point['flow_unit']))
                            pass
                        else:
                            writing.append(dict(set='DionexLC', subset=key, parameter=param, value=hplc[key][param]))

                else:
                    if key in wantDionex or not wantDionex:
                        writing.append(dict(set='DionexLC', subset='main', parameter=key, value=hplc[key]))

            self.logs.log.info('Dionex LC has %i parameters' % len(hplc))
        elif hplc['lctype'] == 'Unknown':
            self.logs.log.info('Unknown or missing HPLC')
            writing.append(dict(set='HPLC', subset='main', parameter='lctype', value='Unknown'))

        # now output the data to HDF5 file
        try:
            writing.sort(key=lambda x: (x['set'], x['subset']))
            self.hdf.appendRows('/rawdata/parameters', writing)

            if profile:
                self.hdf.createTable('rawdata', 'profile', profileTable)
                self.hdf.appendRows('/rawdata/profile', profile)
        except ExHa.HDFmissingDataError, czEx:
            ExHa.addContext(czEx, 'Error in writeLCparameters')
            raise
        except ExHa.HDF5consistancyError, czEx:
            ExHa.addContext(czEx, 'Error in writeLCparameters')
            raise

    # noinspection PyUnusedLocal
    def writeMSparameters(self, ms, wantParam, wantSubsets, order):
        """
        @brief organises the MS parameters for writing to HDF5 file
        @param ms <dictionary>: containing all the MS parameters
        @param wantParam <list>: containing the parameter names to be utilised
        @param wantSubsets <list>: containing the subsets to be utilised
        @param order <list>: containing the order that subsets should be written
        @return:
        """
        hdf = self.hdf
        writing = []

        for subset in order:
            if subset == 'Global Parent Masses' or subset == 'Global Reject Masses':
                # deal with inclusion/reject mass list
                # add data to main parameters table first
                writing.append(dict(set='MS instrument', subset=subset, parameter='numpts', value=len(ms[subset])))
                # mainset, subset, 'numpts', len(params[subset]))

                # only add inclusion table if there is any data
                if ms[subset]:
                    if 'Reject' in subset:
                        hdf.createTable('rawdata', 'reject', 'Reject')
                        reject = []
                        for i in ms[subset]:
                            reject.append(dict(start=i['Start (min)'], end=i['End (min)'], mass=i['Mass']))
                        self.hdf.appendRows('/rawdata/reject', reject)
                    else:
                        hdf.createTable('rawdata', 'inclusion', 'Inclusion')
                        incl = []
                        for i in ms[subset]:
                            d = dict(ms_mass=i['MS Mass'], start=i['Start (min)'], end=i['End (min)'],
                                     faims_cv=i['MS FAIMS CV'], ms_col_energy=i['MS Normalized Collision Energy'],
                                     ms_charge=i['MS Charge State'], ms_intensity=i['MS Intensity Threshold'],
                                     ms2_mass=i['MS2 Mass'], ms2_col_energy=i['MS2 Normalized Collision Energy'],
                                     name=i['Name'])
                            incl.append(d)

                        self.hdf.appendRows('/rawdata/inclusion', incl)
            elif isinstance(ms[subset], dict):
                for param in ms[subset]:
                    writing.append(dict(set='MS instrument', subset=subset, parameter=param, value=ms[subset][param]))
            elif isinstance(ms[subset], list):
                if subset == 'Additional Microscans':
                    writing.append(dict(set='MS instrument', subset='main', parameter=subset, value=str(ms[subset])))
                else:
                    for param in ms[subset]:
                        writing.append(dict(set='MS instrument', subset=subset, parameter=param,
                                            value=ms[subset][param]))
            else:
                if subset in wantParam or not wantParam:
                    writing.append(dict(set='MS instrument', subset='main', parameter=subset, value=ms[subset]))
                pass

        self.hdf.appendRows('/rawdata/parameters', writing)

    # noinspection PyUnusedLocal
    def writeExactiveMSparameters(self, ms, wantParam, wantSubsets, order):
        """
        @brief organises the MS parameters for writing to HDF5 file
        @param ms <dictionary>: containing all the MS parameters
        @param wantParam <list>: containing the parameter names to be utilised
        @param wantSubsets <list>: containing the subsets to be utilised
        @param order <list>: containing the order that subsets should be written
        @return:
        """

        exactive = {
            'INCLUSION LIST': dict(tableName='inclusion', tableClass='ExactiveInclusion',
                                   method=self.writeExactiveInclusion),
            'EXCLUSION LIST': dict(tableName='exclusion', tableClass='ExactiveInclusion',
                                   method='self.writeExactiveInclusion'),
            'NEUTRAL LOSSES': dict(tableName='neutrallosses', tableClass='ExactiveInclusion',
                                   method='self.writeExactiveInclusion'),
            'LOCK MASSES': dict(tableName='lockmass', tableClass='ExactiveLockMass',
                                method=self.writeExactiveLockMass)
        }

        writing = []
        main = [x for x in ms if not isinstance(ms[x], dict)]
        for param in main:
            if param in wantParam or not wantParam:
                writing.append(dict(set='MS instrument', subset='main', parameter=param, value=ms[param]))

        for subset in wantSubsets:
            # only iterate over the wanted subset data
            if subset == 'EXPERIMENT':
                # main experiment data: can be organised in multiple experiments
                exp = ms[subset]
                if len(ms[subset]) == 1:
                    # single experiment, but with greater depth
                    exp = exp[exp.keys()[0]]

                main = [x for x in exp if not isinstance(exp[x], dict)]
                for param in main:
                    writing.append(dict(set='MS instrument', subset=subset, parameter=param, value=exp[param]))

                expSubsets = [x for x in exp if isinstance(exp[x], dict)]
                for expSub in expSubsets:
                    for param in exp[expSub]:
                        writing.append(dict(set='MS instrument', subset=expSub, parameter=param,
                                            value=exp[expSub][param]))

            elif subset in exactive:
                # other subsets have a number of entries
                entries = ms[subset]['entries']
                writing.append(dict(set='MS instrument', subset=subset, parameter='entries', value=entries))
                if entries > 0:
                    # there is data for this subset
                    exactive[subset]['method'](ms[subset], exactive[subset])
            else:
                for param in ms[subset]:
                    writing.append(dict(set='MS instrument', subset=subset, parameter=param, value=ms[subset][param]))

        self.hdf.appendRows('/rawdata/parameters', writing)

    def writeExactiveLockMass(self, lockMasses, tableData):
        writing = []
        for mass in lockMasses['data']:
            writing.append(dict(mass=mass['Mass'], polarity=mass['Polarity'], start=mass['Start'], end=mass['End'],
                                comment=mass['Comment'][:100]))

        self.hdf.createTable('rawdata', tableData['tableName'], tableData['tableClass'])

        self.hdf.appendRows('/rawdata/' + tableData['tableName'], writing)

    def writeExactiveInclusion(self, inclusion, tableData):

        writing = []
        for ion in inclusion['data']:
            tmpDict = dict(ms_mass=float(ion['Mass']), polarity=ion['Polarity'][:3], comment=ion['Comment'][:100])

            if ion['Start'].strip():
                tmpDict['start'] = float(ion['Start'])

            if ion['End'].strip():
                tmpDict['end'] = float(ion['End'])

            if ion['CS'].strip():
                tmpDict['charge'] = float(ion['CS'])

            if ion['NCE'].strip():
                tmpDict['col_energy'] = float(ion['NCE'])

            writing.append(tmpDict.copy())

        self.hdf.createTable('rawdata', tableData['tableName'], tableData['tableClass'])

        self.hdf.appendRows('/rawdata/' + tableData['tableName'], writing)

    def writeTuneParameters(self, tune, wantParam, wantSubset):
        """
        @brief organises the Tune parameters for writing to HDF5 file
        @param tune <dictionary>: containing the tune parameters
        @param wantParam <list>: containing the parameter names to be utilised
        @param wantSubsets <list>: containing the subsets to be utilised
        @return:
        """
        writing = []
        for subset in tune:
            if isinstance(tune[subset], dict):
                if subset in wantSubset:
                    for param in tune[subset]:
                        writing.append(dict(set='Tune', subset=subset, parameter=param, value=tune[subset][param]))
            else:
                if subset in wantParam or not wantParam:
                    writing.append(dict(set='Tune', subset='main', parameter=subset, value=tune[subset]))
                pass
        self.hdf.appendRows('/rawdata/parameters', writing)

    def writeUnitData(self, units, probUnits):
        """
        @brief organises and writes the Unit data to HDF5
        @param units <dictionary>: containing the Unit data
        @param probUnits <integer/dictionary>: 0 if no problem units or dictionary containing the problem units
        @return:
        """

        writing = []
        # only write the MS/MS unit index = 1
        for order, unit in enumerate(units[1]):
            tmp = unit.copy()
            tmp['unit'] = 1
            tmp['order'] = order + 1
            se = str(tmp.pop('scans'))
            tmp['scanevents'] = se[:100]
            writing.append(tmp)

        if probUnits:
            for p in range(1, len(probUnits)):
                for order, unit in enumerate(probUnits[p]):
                    tmp = unit.copy()
                    tmp['unit'] = -p
                    tmp['order'] = order + 1
                    se = str(tmp.pop('scans'))
                    tmp['scanevents'] = se[:100]
                    writing.append(tmp)

        self.hdf.appendRows('/rawdata/units', writing)

    def writeIsotopes(self, isotopes):

        writing = []
        method_id = isotopes['meth_id']
        # batch_id = isotopes['batch_id']
        error = isotopes['match']

        keys = isotopes['quantmasses'].keys()
        keys.sort()
        for key in keys:
            for iso in isotopes['quantmasses'][key]:
                amino = iso.get('amino', '')
                if len(iso['name']) > 10:
                    name = iso['name'][:10]
                    msg = 'Only 10 characters allowed, truncating isotope name <%s> to <%s>' % (iso['name'], name)
                    self.logs.log.warning(msg)
                else:
                    name = iso['name']
                writing.append(dict(iso_id=iso['id'], name=name, mz=iso['mass'], error=error,
                                    intmz=int(iso['mass'] + 0.5), method_id=method_id, amino=amino))

        self.hdf.appendRows('/rawdata/isotopes', writing)

    def loadSpectrum(self, spec, lastsurvey, spRow):
        """
        @brief loads the spectrum data and spectrum related parameters, triggers the processing of the spectrum
        @param spec <integer>: the spectrum number
        @param lastsurvey <integer>: the spectrum number of the last survey spectrum
        @param spRow <integer>: the location for the spectrum header data in self.specs (MS/MS data)
        @return points <integer>: number of data points in the spectrum
        """

        # data references
        cfg = self.cfg
        units = self.units
        xRaw = self.raw

        header = xRaw.getSpectrumHeader(spec)
        event = int(header['extra']['Scan Event'])
        if event in cfg.parameters['general']['skipscanevents']:
            # scan event not wanted so skip all processing
            return -1, spRow
        ionp = self.IonProc
        if self.quanmeth:
            correctionFactors = self.quanmeth['correction']
        else:
            correctionFactors = {}
        data = xRaw.getSpectrumData(spec)
        scanevent = int(header['extra']['Scan Event'])
        meth = {}

        if scanevent > 1:
            if header['filter']['scan'] == 'ms':
                # raise MSdataConsistancyError('Spectrum %d error: Using MS scan filter for scan event = %d' %
                #                              (spec, scanevent))

                # allows the error to be skipped
                return -1, spRow

        elif scanevent == 1 and header['filter']['scan'] != 'ms':
            # allows the error to be skipped
            return -1, spRow

        for ev in range(len(units)):
            for frag in range(len(units[ev])):
                if scanevent in units[ev][frag]['scans']:
                    if ev == 0:
                        meth = 'ms'
                    else:
                        if frag == 0:
                            self.unitnum += 1
                        meth = dict(unitnum=self.unitnum, unit=ev, order=frag + 1, scanevent=scanevent,
                                    act=units[ev][frag]['activation'], use=units[ev][frag]['use'])
                    break
            if meth:
                break
        if meth:
            self.currmeth = meth
        else:
            raise ExHa.MSdataConsistancyError('Spectrum %i error: Unable to find event method for event %i' %
                                              (spec, scanevent))

        if meth == 'ms':
            # handle MS survey spectra
            smdata = SpectrumManager(spec, data, header, cfg, self.specproc)
            smdata.process()
            self.hdf.appendRows('/rawdata/spectra', [dict(spec_id=spec, rt=header['rt'],
                                                          type=header['filter']['scan'])])

            self.writeSpectrumParameters(spec, header, cfg.parameters['wanted']['extra_param'])

            # save spectrum data to generate XIC information
            self.writeXICbinData(spec, header['rt'], smdata.ions)

            self.mssurvey[spec] = dict(rt=header['rt'])

            return smdata, spRow
        else:
            # check if the header has an appropriate monoisotopic mz value
            mono = float(header['extra']['Monoisotopic M/Z'])
            recalc = 0
            if mono == 0:
                # no value set so should look in MS spectrum
                self.logs.log.debug('Monoisotopic value is zero.')
                recalc = 1
            elif mono > float(header['filter']['setmass1']) + 0.5:
                # probabaly a miss calculated data value
                self.logs.log.debug('Monoisotopic value too high.')
                recalc = 1
            elif header['extra']['Charge State'] == '0':
                # probabaly a miss calculated data value
                self.logs.log.debug('Charge state = 0')

            if recalc:
                prec = self.surveydata.findPrecursor(float(header['filter']['setmass1']))
                if prec == 0:
                    msg = 'Setmass (%s) not found in Survey (%s): using setmass with charge state 2.'
                    self.logs.log.debug(msg % (header['filter']['setmass1'], self.surveydata.specID))
                    # no matching precursor could be found use setmass and charge = 2+
                    header['extra']['Monoisotopic M/Z'] = header['filter']['setmass1']
                    header['extra']['Charge State'] = '2'
                else:
                    # precursor was found
                    header['extra']['Monoisotopic M/Z'] = '%.4f' % prec[0]
                    header['extra']['Charge State'] = str(prec[1])

            # header['prec_mz'] = self.surveydata.findion(float(header['extra']['Monoisotopic M/Z']))

            # handle MS/MS spectra
            smdata = SpectrumManager(spec, data, header, cfg, self.specproc, lastsurvey)
            # usehcd = 0

            # do the quantitation
            sumrepint = 0
            sumreparea = 0

            ionFillTime = float(header['extra']['Ion Injection Time (ms)'])
            if self.quan and 'Q' in meth['use']:

                quantIons = smdata.findReporterIons(self.isotopes, self.lastQuant, self.coalescence)

                if quantIons:
                    # do isotopecorrection at this point - tobmat 20141104
                    # we make a dictionary from the quanions array: this is then supplied to the
                    # doIsotoepCorrection function
                    rawQuantVals = {}

                    for labelID in quantIons['isos']:
                        ion = quantIons['isos'][labelID]
                        ion['area'] = ion['inten'] * ionFillTime
                    ionp.doIsotopeCorrection(correctionFactors, quantIons['isos'])
                    sumrepint = quantIons['sumint']
                    sumreparea = quantIons['sumint'] * ionFillTime
                    self.writeQuanData(spec, quantIons, lastsurvey, ionFillTime, self.quan)
                    self.offsets.append(quantIons['offset_mean'])
                    if quantIons['num'] >= self.lastQuant[0]:
                        self.lastQuant = (quantIons['num'], quantIons['offset_mean'])

            # create base entry for the header array - no id/quan spec info
            ns = (spec, meth['unitnum'], meth['order'], meth['scanevent'], float(header['extra']['Monoisotopic M/Z']),
                  0, 0, int(header['extra']['Charge State']), 0, float(header['filter']['setmass1']), lastsurvey, 0, 0,
                  header['rt'], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  header['filter']['frag1'], float(header['filter']['energy1']), sumrepint, sumreparea)

            spRow += 1

            self.specs.put([spRow], ns)
            self.hdf.appendRows('/rawdata/spectra', [dict(spec_id=spec, rt=header['rt'],
                                                          type=header['filter']['scan'])])

            if smdata.ions:
                self.writeMSMSions(spec, smdata.ions, header['rt'])
            self.writeSpectrumParameters(spec, header, cfg.parameters['wanted']['extra_param'])

            return 0, spRow

    def writeSpectrumParameters(self, specID, parameters, wanted):
        """
        @brief organises and sends the spectrum parameter data for adding to the HDF5 file
        @param specID <integer>: the identifier for the spectrum
        @param parameters <dictionary>: containing the parameter data for one spectrum
        @param wanted <list>: containing the names of the 'extra' parameters that are wanted
        @return:
        """

        writing = []

        # first add selected extra parameters
        for want in wanted:
            writing.append(dict(spec_id=specID, set='extra', parameter=want, value=parameters['extra'].get(want, '')))

        # then add all the filter parameters
        for key in parameters['filter']:
            value = parameters['filter'][key]
            if value:
                writing.append(dict(spec_id=specID, set='filter', parameter=key, value=value))

        self.hdf.appendRows('/rawdata/specparams', writing)

    def writeMSMSions(self, specID, data, rt):
        """
        @brief organises and sends the MS/MS spectrum ions for adding to the HDF5 file
        @param specID <integer>: the identifier for the spectrum
        @param data <list>: contains the ion data, if the data has been processed then the ion data is a dictionary
        otherwise a tuple
        @return:
        """

        writing = []
        if isinstance(data[0], dict):
            # this is processed MS/MS data
            for ion in data:
                writing.append(dict(spec_id=specID, mz=ion['mz'], mzraw=ion['mz'], inten=ion['topint'], rt=rt))
        else:
            # non-processed MS/MS data - collected in stick mode
            for ion in data:
                writing.append(dict(spec_id=specID, mz=ion[0], inten=ion[1], rt=rt))

        self.hdf.appendRows('/rawdata/ions', writing)

    def writeXICbinData(self, specID, rt, iondata):
        """
        @brief organises and sends the MS1 centroided ions for adding to the HDF5 file
        @param specID <integer>: the identifier for the spectrum
        @param rt <float>: the retention time of the spectrum
        @param iondata <ndarray>: containing the centroided peaks from MS1 data
        @return:
        """

        writing = []
        for ion in iondata:
            writing.append(dict(bin=ion['bin'], specid=specID, rt=rt, mz=ion['mz'], mzraw=ion['mz'],
                                inten=ion['inten'], area=ion['area']))

        self.hdf.appendRows('/rawdata/xicbins', writing)

    def writeNoise(self, specID, noiseData):
        """
        @brief organises and sends the spectrum noise data for adding to the HDF5 file
        @param specID <integer>: the identifier for the spectrum
        @param noiseData <list>: containing noise data in a tuple of (mz, noise, baseline)
        @return:
        """

        writing = []
        for n in noiseData:
            writing.append(dict(spec_id=specID, mz=n[0], noise=n[1], baseline=n[2]))

        self.hdf.appendRows('/rawdata/noise', writing)

    def writeQuanData(self, specID, quanions, survey, ionFillTime, repions):
        """
        @brief organises and sends the spectrum noise data for adding to the HDF5 file
        @param specID <integer>: the identifier for the spectrum
        @param quanions <dictionary>: of all possible quantification ions and their detected ions
        @param survey <integer>: identifier for the survey spectrum prior to the MS/MS event
        @param ionfill <float>: the time taken to fill the trap
        @param repions <dictionary>: theoretical data for the quantification isotopes
        @return:
        """

        writing = []
        isos = quanions['isos']
        orderedIDs = isos.keys()
        orderedIDs.sort()
        for labelID in orderedIDs:
            if isos[labelID]['mz'] > 0:
                ion = isos[labelID]
                writing.append(dict(spec_id=specID, isolabel_id=labelID, survey_id=survey, mzdiff=ion['delta'],
                                    area=ion['inten'] * ionFillTime, inten=ion['inten'], coalescence=ion['coalescence'],
                                    ppm=ion['delta'] / repions['quantmasses'][labelID]['mass'] * 1e6,
                                    corrected=ion['corrected']))
        if writing:
            self.hdf.appendRows('/rawdata/quan', writing)

    def createIndexes(self):
        """
        @brief organises the creation of indexes on appropriate tables in the hdf5 file
        @return:
        """

        hdf = self.hdf
        hdf.indexTable('/rawdata/ions', ['spec_id', 'mz', 'inten'])
        hdf.indexTable('/rawdata/spectra', ['spec_id'])
        hdf.indexTable('/rawdata/specparams', ['spec_id'])
        hdf.indexTable('/rawdata/xicbins', ['specid', 'bin', 'rt', 'inten'])
        hdf.indexTable('/rawdata/parameters', ['set'])
        hdf.indexTable('/rawdata/noise', ['spec_id'])
        pass

    # @staticmethod
    # def getIsootopeData(peak, theory):
    #     # theoretical data
    #     expect = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0}
    #     for i in range(min(4, len(theory))):
    #         expect[i] = theory[i]
    #     # found data
    #     found = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0}
    #     for isodata in peak['isotopes']:
    #         iso = isodata[0] - peak['offset']
    #         if iso in [0, 1, 2, 3]:
    #             found[iso] = isodata[1]
    #
    #     return found, expect

    def processXICs(self):
        """
        @brief extracts the XIC data for the MS/MS events or the provided mass-time list
        @param watch <object>: stopwatch object that records timings for the analysis
        """

        # methods
        global isomodel
        deisotopexicdata = self.deisotopeXICdata
        calcS2I = self.calcS2I
        selpeak = self.selBestFittingPeak
        writeheader = self.writeHeader

        xic = XIC(self.hdf, self.cfg)

        # references to data
        specs = self.specs
        cfg = self.cfg
        ppm = cfg.parameters['general']['tolppm']
        mda = cfg.parameters['general']['tolmda']
        wfid = self.wf_id

        # open shelve object
        precursorShelf = shelve.open(self.tempFile, protocol=2)
        ordered = []
        for key in precursorShelf.iterkeys():
            ordered.append((int(key), key))

        ordered.sort()

        self.logs.log.info('Processing XIC data for %d units' % len(ordered))

        # index the spectra
        specindex = {}
        for i, s in enumerate(specs):
            specindex[s['spec']] = i

        pBar = progBar.ProgressBar(widgets=progBar.name_widgets, maxval=len(ordered), name='XIC analysis').start()
        # cnt = Counter(total=len(ordered))

        cnt = 0
        try:
            for unitkey in ordered:
                sw = Stopwatch()

                cnt += 1
                posspeaks = []
                pBar.update(cnt)

                unit = precursorShelf[unitkey[1]]

                if unit['charge'] == 0:
                    unit['precmz'] = unit['precmz_surv']
                    unit['offset'] = unit['offset_surv']
                    unit['precmz_xic'] = 0.0
                    unit['offset_xic'] = 0
                    for sp in unit['msms']['specs']:
                        unit['msms']['specs'][sp]['s2i'] = 0.0
                        unit['msms']['specs'][sp]['c12'] = 0.0
                    unit['early']['s2i'] = 0.0
                    unit['late']['s2i'] = 0.0
                    unit['early']['c12'] = dict(fit=0, base=0, c12=0, inten=0, mz=0, norm=0, offset=0, smooth=0)
                    unit['late']['c12'] = dict(fit=0, base=0, c12=0, inten=0, mz=0, norm=0, offset=0, smooth=0)
                    continue

                # test for the presence of the bin range from the precursor
                unitmz = unit['basemass']

                mzerror = max((mda, unitmz * ppm))

                # find the isotope ratios for this ion
                molmass = unitmz * unit['charge']
                isoset = isomodel.findIsotopeSet(molmass, 2)

                # find range for isotpeblock generation
                lowmass = unit['msms']['isoions'][0][1]
                maxiso = abs(unit['msms']['isoions'][0][0])
                for u in unit['msms']['isoions']:
                    if u[0] == maxiso:
                        highmass = u[1]
                        break

                xic.setXICblock(lowmass, highmass, unit['xicmin'], unit['xicmax'])
                # extract XIC data and extract peaks for precursor ion
                xicLen = xic.createXIC(unitmz, mzerror)

                if xicLen:
                    posspeaks = xic.findPeaksInXIC(unit['xicmin'], unit['xicmax'], unit['baseoffset'])
                    if posspeaks:
                        # now look for isotope data

                        # loop through all other isos skipping the precursor
                        for iso in unit['msms']['isoions']:
                            # stop when maxiso is reached, minimising the loading of unrequired data
                            if iso[0] > maxiso:
                                break
                            if iso[0] == unit['baseoffset']:
                                # skip the precursor
                                continue
                            isomz = iso[1]

                            xicLen = xic.createXIC(isomz, mzerror)
                            if xicLen:
                                newpeaks = xic.findPeaksInXIC(unit['xicmin'], unit['xicmax'], iso[0])

                                if newpeaks:
                                    xic.linkOverlappingPeaks(posspeaks, newpeaks, iso)

                deiso = []
                for peak in posspeaks:
                    peak.sort(key=lambda x: x['offset'])
                    nextpeak = deisotopexicdata(peak, isoset['inten'], unit['charge'])
                    if nextpeak:
                        if nextpeak['fit'] < unit['surv_fit'] or nextpeak['offset'] == unit['offset_surv']:
                            nextpeak['deltart'] = nextpeak['rt'] - unit['rt']
                            deiso.append(nextpeak.copy())

                # order the peaks by closeness to MS/MS event
                deiso.sort(key=lambda x: abs(x['deltart']))

                match = selpeak(deiso)
                # sw.rec('selection')
                if match:
                    # found a isotope match
                    unit['xicpeak'] = match.copy()
                    unit['precmz_xic'] = match['peakmz']
                    unit['offset_xic'] = match['offset']
                    unit['isos'] = isoset['inten']

                else:
                    unit['precmz_xic'] = 0.0
                    unit['offset_xic'] = -100

                # now find the best peak to match to the MS/MS event
                calcS2I(unit)

                # set the precursor m/z
                if unit['precmz_xic'] == 0:
                    # no xic data use survey
                    unit['precmz'] = unit['precmz_surv']
                    unit['offset'] = unit['offset_surv']
                    unit['type'] = 'survey'
                elif unit['precmz_surv'] == 0:
                    # no survey data use xic
                    unit['precmz'] = unit['precmz_xic']
                    unit['offset'] = unit['offset_xic']
                    unit['type'] = unit['xicpeak']['type']
                elif unit['offset_xic'] == unit['offset_surv']:
                    # survey and xic have the same offset so use the XIC mass
                    unit['precmz'] = unit['precmz_xic']
                    unit['offset'] = unit['offset_xic']
                    unit['type'] = unit['xicpeak']['type']
                else:
                    # use the tightest fitting data
                    if unit['surv_fit'] < unit['xicpeak']['fit']:
                        # survey is better
                        unit['precmz'] = unit['precmz_surv']
                        unit['offset'] = unit['offset_surv']
                        unit['type'] = 'survey'
                    else:
                        # XIC is better
                        unit['precmz'] = unit['precmz_xic']
                        unit['offset'] = unit['offset_xic']
                        unit['type'] = unit['xicpeak']['type']

                unitSpecs = unit['msms']['specs'].keys()
                unitSpecs.sort()

                if unit['precmz'] == 0:
                    # still no match: use monoisotopicMZ
                    unit['precmz'] = specs[specindex[unitSpecs[0]]]['monomz']
                    unit['type'] = 'monomz'

                for spec_id in unitSpecs:
                    writeheader(unit, specs[specindex[spec_id]])
                # sw.rec('writing')
                sw.stop()
                self.logs.log.debug('%i\t%.5f\t%.2f\t%s' % (unitkey[0], unit['basemass'], unit['rt'],
                                                            sw.result[-1]['split']))
            pBar.finish()

        except:
            # other exceptions
            trace = repr(traceback.format_exception(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
            return {'code': 1, 'error': 'pyMSsafe XIC processing error (unit %s): %s' % (unitkey[1], trace)}

        return {'code': 0}

    def writeHeader(self, unit, msms):

        header = dict(spec_id=msms['spec'], unit_id=msms['unit'], order=msms['order'], scanevent=msms['scanevent'],
                      precmz=unit['precmz'], precmzraw=unit['precmz'], precmz_xic=unit['precmz_xic'],
                      monomz=msms['monomz'], precmz_surv=msms['precmz_surv'],
                      charge=msms['preccharge'], precsd=msms['precsd'], setmass=msms['setmass'], rt=msms['rt'],
                      survey_spec=msms['survey'], quan_spec=msms['quanspec'], id_spec=msms['idspec'],
                      fragmeth=msms['fragmeth'].upper(), fragenergy=msms['fragenergy'],
                      sumrepint=msms['sumrepint'], sumreparea=msms['sumreparea'],
                      s2i=unit['msms']['specs'][msms['spec']]['s2i'], s2iearly=unit['early']['s2i'],
                      s2ilate=unit['late']['s2i'], thresh=unit['msms']['specs'][msms['spec']]['thresh'],
                      threshearly=unit['early']['thresh'], threshlate=unit['late']['thresh'],
                      c12=unit['msms']['specs'][msms['spec']]['c12'], c12early=unit['early']['c12']['inten'],
                      c12late=unit['late']['c12']['inten'],
                      peaktype=unit['type'], precinten=unit['msms']['specs'][msms['spec']]['c12'])

        if unit['precmz_xic'] > 0:
            extraDic = dict(scanapex=unit['xicpeak']['specid'], rtapex=unit['xicpeak']['rt'],
                            area=unit['xicpeak']['area'], inten=unit['xicpeak']['inten'],
                            smooth=unit['xicpeak']['smooth'], integ=unit['xicpeak']['integ'],
                            fwhm=unit['xicpeak']['rt50right'] - unit['xicpeak']['rt50left'],
                            c13iso=unit['xicpeak']['c13'], precinten=unit['xicpeak']['inten'])

            header.update(extraDic)

        self.hdf.appendRows('/rawdata/msmsheader', [header])

    @staticmethod
    def selBestFittingPeak(deiso):
        if not deiso:
            return 0
        deiso.sort(key=lambda x: x['fit'])

        for pk in deiso:
            if pk['c13ok'] > 0:
                match = pk.copy()
                match['type'] = 'XICclust'
                return match
        match = deiso[0].copy()
        match['type'] = 'XICsingle'
        return match

    def deisotopeXICdata(self, xicdata, isoratios, charge):
        """
        @brief creates the best isotope data from a list of XIC peaks
        """

        isotopes = []
        ions = {}
        for xic in xicdata:
            ions[xic['offset']] = (xic['peakmz'], xic['smooth'], xic['inten'], xic['offset'])
            isotopes.append((xic['offset'], xic['smooth'], xic['peakmz']))

        fit = self.isoFitLeastSQ(isoratios, ions, charge, 1)
        for x in xicdata:
            if x['offset'] == fit[3]:
                b = x.copy()
                b['c13'] = fit[1].count(1) - 1
                if fit[1][1] == 0:
                    b['c13ok'] = 0
                elif fit[1][1] == 1:
                    if b['c13'] > 1:
                        b['c13ok'] = 2
                    else:
                        b['c13ok'] = 1
                b['isotopes'] = isotopes
                b['norm'] = fit[2]
                b['fit'] = fit[0]
                return b

    def updateMSMS(self, start, stop):
        """
        @brief updates the id data for the specified range of spectra
        @param start <integer>: low end of the range to update
        @param stop <integer>: high end of the range to update
        """
        ishcd = self.ishcd
        # find the data that needs updating
        new = {}
        specs = self.specs
        for spec in range(start, stop):
            if specs[spec]['monomz'] > 0:
                new.setdefault(specs[spec]['setmass'], []).append(spec)

        # now test the data
        for mz in new:
            if ishcd == 1 or len(new[mz]) == 1:
                for sp in new[mz]:
                    specs[sp]['quanspec'] = specs[sp]['spec']
                    specs[sp]['idspec'] = specs[sp]['spec']
            elif len(new[mz]) % ishcd == 0:
                for j in range(0, len(new[mz]), 2):
                    first = new[mz][j]
                    second = new[mz][j + 1]
                    if specs[first]['quanspec'] == 0:
                        specs[first]['quanspec'] = specs[second]['quanspec']
                        specs[second]['idspec'] = specs[first]['idspec']
                    elif specs[second]['quanspec'] == 0:
                        specs[second]['quanspec'] = specs[first]['quanspec']
                        specs[first]['idspec'] = specs[second]['idspec']
            else:
                lspecs = []
                for j in new[mz]:
                    lspecs.append(specs[spec]['spec'])

                msg = 'Unexpected MS/MS pattern %i events together not %i (%s)' % (len(new[mz]), max(1, ishcd),
                                                                                   ','.join(lspecs))
                raise ExHa.MSdataConsistancyError(msg)

    def calcBackgroundIsotopes(self, specID, survey):
        """
        @brief sets the s2i data from the survey data for the MS/MS spectrum given by id
        @param id <integer>: MS/MS spectrum id
        @param survey <object>: SpectrumManager object containing MS survey data
        @param setby <function>: method by which the "best" signal is selected
        @return background <float>: sum of all the background ions
        @return isoions <list>: containing the data for potential isotope ions
        """
        # data references
        cfgS2I = self.cfg.parameters['sig2interf']
        winfull = cfgS2I['winfull']
        winscaled = cfgS2I['winscaled']
        neutron = self.cfg.parameters['general']['neutron']
        maxIso = self.cfg.parameters['deisotoping']['max_isotopes']
        intenScale = cfgS2I['intenscale']
        msions = survey.ions

        msms = self.specs[specID]
        mz = float(msms['setmass'])
        isomz = float(msms['monomz'])
        low = mz - winscaled
        high = mz + winscaled
        lowfull = mz - winfull
        highfull = mz + winfull
        # isorange = winfull * 2
        mda = self.cfg.parameters['general']['tolmda']
        ppm = self.cfg.parameters['general']['tolppm']
        error = max(mda, isomz * ppm)
        matchby = cfgS2I['matchby']

        # extract ions from the isolation window
        ions = [x for x in msions if low <= x['mz'] <= high]

        background = 0
        isoions = {}

        # calculate series of ions corresponding to neutron spacing
        # use a range from -max_offset to +maxisotopes
        isoset = isomodel.findIsotopeSet(mz * msms['preccharge'], 2)
        isomzs = [(i, isomz + i * neutron / msms['preccharge']) for i in range(1 - len(isoset['inten']), maxIso)]

        if ions:
            # calculate the background
            iso = 0
            stop = 0
            for ion in ions:
                if ion['mz'] < lowfull or ion['mz'] > highfull:
                    inten = ion['inten'] * intenScale
                else:
                    inten = ion['inten']
                background += inten

                # skip to the next isotope if needed
                while ion['mz'] > isomzs[iso][1] + error:
                    iso += 1
                    if iso >= len(isomzs):
                        stop = 1
                        break
                if stop:
                    break
                if ion['mz'] < isomzs[iso][1] - error:
                    # passed
                    pass
                elif isomzs[iso][0] in isoions:
                    # has ion already present for that isotope
                    if matchby == 'mass':
                        # use the closest in mass
                        if abs(ion['mz'] - isomzs[iso][1]) < abs(isoions[isomzs[iso][0]][0] - isomzs[iso][1]):
                            isoions[isomzs[iso][0]] = (ion['mz'], ion['inten'], inten, isomzs[iso][0])
                    else:
                        if isoions[isomzs[iso][0]][1] < ion['inten']:
                            isoions[isomzs[iso][0]] = (ion['mz'], ion['inten'], inten, isomzs[iso][0])
                else:
                    isoions[isomzs[iso][0]] = (ion['mz'], ion['inten'], inten, isomzs[iso][0])

        return background, isoions, isomzs

    def calcS2I(self, unit):
        """
        @brief calculates the s2i for the ions
        @param unit <dictionary>: cotaining the pre and post MS/MS survey data
        @param cfg <object>: contains configuration information
        """
        # methods
        interpolate = self.interpolate
        # for each spectrum in the unit
        for sp in unit['msms']['specs']:
            spec = unit['msms']['specs'][sp]

            # calculate the s2i from the early and late spectra
            for pos in ['early', 'late']:
                c12 = (0, 0, 0, 0)
                if unit[pos]['bkgd'] == 0:
                    # if background is 0 then no signal everything = 0
                    unit[pos]['s2i'] = 0.0
                else:
                    signal = 0.0
                    # check if isotopes found
                    if 'foundisos' in unit[pos]:
                        # find the 122C ion
                        offset = unit.get('offset', 0)
                        if offset in unit[pos]['foundisos']:
                            c12 = unit[pos]['foundisos'][offset]

                        isos = unit[pos]['foundisos'].keys()
                        isos.sort()
                        last = offset
                        for iso in isos:
                            if iso < offset:
                                continue
                            elif iso <= last + 1:  # use only adjacent ions
                                signal += unit[pos]['foundisos'][iso][2]
                                last = iso

                    else:
                        pass

                    unit[pos]['s2i'] = signal / unit[pos]['bkgd']

                if pos == 'early':
                    earlyc12 = c12[2]
                else:
                    latec12 = c12[2]

                if unit[pos]['s2i'] > 1:
                    pass
            spec['s2i'] = interpolate(unit['early']['s2i'], unit['early']['rt'],
                                      unit['late']['s2i'], unit['late']['rt'], spec['rt'])
            spec['c12'] = interpolate(earlyc12, unit['early']['rt'],
                                      latec12, unit['late']['rt'], spec['rt'])
            pass


class SpectrumManager():
    def __init__(self, specID, data, header, cfg, specproc, lastsurvey=0):
        """
        @param id <integer>: the spectrum number
        @param data <list>: containing the ion tuples (m/z, inten) for the entire spectrum
        @param header <dictionary>: containing all the spectrum specific parameters
        @param cfg <config object>: containing all the configuration parameters
        @param specproc <SpectrumProcessor Object>: used to process the data
        @param lastsurvey <integer>: the id of the last survey (only used for MS/MS data)
        """
        self.specID = specID
        self.header = header
        self.cfg = cfg
        self.type = header['filter']['scan']
        self.noise = None
        self.processor = specproc
        self.binsPerDa = cfg.parameters['xic']['binsperda']
        self.mininten = 0

        if header['filter']['scan'] in ('ms2', 'ms3'):
            self.survey = lastsurvey
            self.data = 0
            if header['filter']['centroid'] == 'c':
                self.ions = data
                self.isProfile = False
            elif header['filter']['centroid'] == 'p':
                self.data = 0
                self.ions, minint = self.processor.process(data, pick='gaussian', isMS2=1)
                self.isProfile = True
        else:
            self.data = data

    def process(self):
        """
        @brief processes the spectrum to generate ions
        @param specproc <object>: containing all the processing methods
        """
        prec = self.binsPerDa
        iondics, self.mininten = self.processor.process(self.data)

        # convert dictionaries to arrays
        iontuples = []
        for ion in iondics:
            iontuples.append((int(ion['mz'] * prec), ion['topmz'], ion['mz'], ion['topint'],
                              ion['area'], ion['logint'], ion['logarea']))

        columns = [('bin', int), ('peaktop', int), ('mz', float), ('inten', float), ('area', float),
                   ('logint', float), ('logarea', float)]
        self.ions = np.array(iontuples, dtype=columns)

    def findReporterIons(self, repIons, lastQuant, coalescence):

        locateTol = repIons['search']
        matchTol = repIons['match']
        flank = max(locateTol * 2, 1.0)
        isotopes = repIons['quantmasses']
        isotopeIDlist = isotopes.keys()
        isotopeIDlist.sort()

        spectrum = self.ions

        if self.isProfile:
            # spectrum processed by pyMsafe: list of dictionaries - set keys
            mass = 'mz'
            inten = 'topint'
        else:
            # spectrum processed by Xcalibur: list of tuples - set indexes
            mass = 0
            inten = 1

        # find the mass range of the reporter ion region
        lowMass = isotopes[isotopeIDlist[0]]['mass']
        highMass = isotopes[isotopeIDlist[0]]['mass']
        for isotopeID in isotopeIDlist:
            if isotopes[isotopeID]['mass'] > highMass:
                highMass = isotopes[isotopeID]['mass']
            elif isotopes[isotopeID]['mass'] < lowMass:
                lowMass = isotopes[isotopeID]['mass']

        # extend by the flanking region
        highMass += flank
        lowMass -= flank

        # extract reporter ion region
        i = 0
        ionTuples = []
        numions = len(spectrum)
        while i < numions and spectrum[i][mass] < highMass:
            if self.isProfile:
                # spectrum processed by pyMsafe: list of dictionaries - set keys
                ionTuples.append((spectrum[i]['mz'], spectrum[i]['topint'],
                                  spectrum[i]['minleft'], spectrum[i]['minright']))
            else:
                # spectrum processed by Xcalibur: list of tuples - set indexes
                ionTuples.append((spectrum[i][0], spectrum[i][1], 0, 0))
            i += 1

        spec = np.array(ionTuples, dtype=[('mz', float), ('inten', float), ('minleft', float), ('minright', float)])
        spec = spec[(spec['mz'] > lowMass)]

        # measure all possible starting points
        deltas = []
        for isotopeID in isotopeIDlist:
            deltas += self.findSeedDeltas(spec, isotopes[isotopeID]['mass'], locateTol)

        # sort the starting deltas and create the startDeltas list of dictionaries
        deltas.sort()
        startDeltas = []
        for d in deltas:
            startDeltas.append(dict(offset_seed=d, isos={}))

        for result in startDeltas:
            offset = result['offset_seed']
            best = {}

            last_id = 0
            errors = []
            sumInt = 0.0
            for iso_id in isotopeIDlist:
                isotope = isotopes[iso_id]
                found = spec[(spec['mz'] >= isotope['mass'] + offset - matchTol) &
                             (spec['mz'] <= isotope['mass'] + offset + matchTol)]

                if len(found) > 0:
                    maxInt = found['inten'].max()
                    for f in found:
                        if f['inten'] == maxInt:
                            best[iso_id] = self.createIsoIonDict(iso_id, isotope['mass'], f)
                            errors.append(best[iso_id]['delta'])
                            sumInt += best[iso_id]['inten']

            best['offset_mean'] = stats.mean(errors)
            best['offset_min'] = min(errors)
            best['offset_max'] = max(errors)
            best['offset_range'] = best['offset_max'] - best['offset_min']
            best['sumint'] = sumInt
            best['num'] = len(errors)
            best['diff2last'] = abs(best['offset_mean'] - lastQuant[1])

            for key in best.keys():
                if key in isotopeIDlist:
                    result['isos'][key] = best[key].copy()
                else:
                    result[key] = best[key]

        startDeltas = sorted(startDeltas, key=lambda x: (-x['num'], -x['sumint'], x['diff2last']))

        # filter the startDeltas to the best results for each number of Reporter ions
        comparison = []
        numbers = []
        if startDeltas:
            maxNum = startDeltas[0]['num']
        else:
            maxNum = 0

        for sd in startDeltas:
            if sd['num'] in numbers or sd['num'] < maxNum - 2:
                continue
            else:
                comparison.append(sd)
                numbers.append(sd['num'])

        comparison = sorted(comparison, key=lambda x: -x['sumint'])
        if comparison:
            bestDelta = comparison[0]

            self.findCoalescence(bestDelta, coalescence, lastQuant)
            if bestDelta['num'] == maxNum:
                betterInt = 0
            else:
                betterInt = 1

            return bestDelta
        else:
            return {}

    def findCoalescence(self, isotopes, coalescence, lastQuant):
        offset = isotopes['offset_mean']
        tol = self.cfg.parameters['quantitation']['coalescence_tol']
        remThresh = self.cfg.parameters['quantitation']['ionremthresh']

        if isotopes['offset_range'] < 2 * tol:
            return

        for pair in coalescence:
            if pair[0] in isotopes['isos']:
                lowIso = isotopes['isos'][pair[0]]
                lowDelta = lowIso['delta'] - offset
            else:
                lowIso = 0
                lowDelta = 0

            if pair[1] in isotopes['isos']:
                highIso = isotopes['isos'][pair[1]]
                highDelta = offset - highIso['delta']
            else:
                highIso = 0
                highDelta = 0

            if lowDelta > tol and highDelta > tol:
                # both ions have moved together
                if lowIso['minright'] == 0:
                    lowIso['coalescence'] = -1
                else:
                    lowIso['coalescence'] = lowIso['minright'] / lowIso['inten']
                if highIso['minleft'] == 0:
                    highIso['coalescence'] = -1
                else:
                    highIso['coalescence'] = highIso['minleft'] / highIso['inten']
            elif lowDelta > tol:
                # lowIso has moved and highIso missing or not moved
                if pair[2] + offset - lowIso['mz'] < lowDelta:
                    # lowIon is fully coalesced
                    lowIso['coalescence'] = 1
                    if highIso and highIso['inten'] / lowIso['inten'] < remThresh:
                        # highIso present and v low intensity - so eliminate
                        del isotopes['isos'][pair[1]]
                        repIonStats = self.calcRepIonStats(isotopes['isos'], lastQuant)
                        isotopes.update(repIonStats)

            elif highDelta > tol:
                # highIso has moved and lowIso missing or not moved
                if highIso['mz'] - pair[2] - offset < highDelta:
                    highIso['coalescence'] = 1
                    if lowIso and lowIso['inten'] / highIso['inten'] < remThresh:
                        # lowIso present and v low intensity - so eliminate
                        del isotopes['isos'][pair[0]]
                        repIonStats = self.calcRepIonStats(isotopes['isos'], lastQuant)
                        isotopes.update(repIonStats)

    @staticmethod
    def calcRepIonStats(repIonSet, lastQuant):
        errors = []
        sumInt = 0.0
        for i in repIonSet:
            errors.append(repIonSet[i]['delta'])
            sumInt += repIonSet[i]['inten']

        minOff = min(errors)
        maxOff = max(errors)
        meanOff = stats.mean(errors)
        statDict = dict(offset_mean=meanOff, offset_min=minOff, offset_max=maxOff, offset_range=maxOff - minOff,
                        sumint=sumInt, num=len(errors), diff2last=abs(meanOff - lastQuant[1]))
        return statDict

    @staticmethod
    def findSeedDeltas(spec, mass, tol):
        deltas = []
        found = spec[(spec['mz'] >= mass - tol) & (spec['mz'] <= mass + tol)]
        for f in found:
            deltas.append(f['mz'] - mass)
            # deltas.append((f['inten'], f['mz'] - mass))
        return deltas

    @staticmethod
    def createIsoIonDict(iso_id, mass, ionData):
        return dict(iso_id=iso_id, mz=ionData['mz'], inten=ionData['inten'], delta=ionData['mz'] - mass,
                    minleft=ionData['minleft'], minright=ionData['minright'], coalescence=0)

    def findPrecursor(self, setmass):
        """
        @brief finds the monoisotopic ion closest to the setmass
        @param setmass <float>: mass used to isolate the precursor for MS/MS
        """
        cfg = self.cfg
        global isomodel
        neutron = cfg.parameters['general']['neutron']
        loghightol = math.log10(cfg.parameters['deisotoping']['high_int_tol'])
        ions = self.ions
        tid = (ions['mz'] > setmass - 5) & (ions['mz'] < setmass + 5)
        local = ions[tid]
        deltas = []
        mindelta = (5, 0)
        for i in local:
            d = i['mz'] - setmass
            if abs(d) < abs(mindelta[0]):
                mindelta = (d, len(deltas))
            deltas.append(d)
        charge = {}
        for chg in range(1, 6):
            isos = {}
            liso = []
            use = 0
            for idx, delta in enumerate(deltas):
                dev = delta - mindelta[0]
                neuts = int(dev / neutron * chg + 0.5)
                ppm = (dev - neuts * neutron / chg) / local[mindelta[1]]['mz']
                if abs(ppm) < 0.000060:
                    isos[neuts] = idx
                    liso.append(neuts)

            if -1 in isos or 1 in isos:
                # have a direct link to precusor
                order = [liso[0]]
                for i in liso:
                    if order[-1] == i - 1:
                        order.append(i)
                    elif 0 not in order:
                        order = [i]
            else:
                order = []

            # extract ions for the series
            if 0 in order:
                charge[chg] = []
                for o in order:
                    charge[chg].append(isos[o])
        # now check each series for their isotope patterns
        prec = []
        for chg in charge.keys():
            found = charge[chg]
            isoratios = isomodel.findIsotopeSet(setmass * chg)
            ok = [found[0]]
            c12isoint = local[found[0]]['inten']

            for i in range(1, min(len(found), len(isoratios['inten']) - isomodel.c12)):
                change = math.log10(local[found[i]]['inten'] / c12isoint / isoratios['inten'][isomodel.c12 + i])
                if abs(change) <= loghightol:
                    # intensities OK
                    ok.append(found[i])
                else:
                    ok.append(-1)
                    break

            base = False
            for idx in range(len(ok)):
                if found[idx] == mindelta[1] and ok[idx] != -1:
                    base = True
                elif ok[idx] == -1:
                    break

            if idx >= 2 and base:
                prec.append((local[found[0]]['mz'], chg))

        if len(prec) == 1:
            return prec[0]
        elif len(prec) > 1:
            return prec[1]
        else:
            return 0

    def getThreshold(self, ion):

        noise = self.noise
        p = 0
        try:
            while noise[p][0] <= ion:
                p += 1
        except IndexError:
            return noise[-1][1]

        if p == 0:
            return noise[p][1]
        else:
            high = noise[p]
            low = noise[p - 1]
            thresh = ((high[1] - low[1]) * (ion - low[0]) / (high[0] - low[0])) + low[1]
            return thresh
