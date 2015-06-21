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

This file controls the output of MS/MS data into Mascot readable .mgf files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

# python modules
import numpy as np
import sys
from pathlib import Path

# CommonUtils modules

sys.path.insert(0, '..')
from CommonUtils.ConfigManager import pyMSsafeConfigManager
from CommonUtils.tools import *
from CommonUtils.hdf5Base import hdf5Base
import CommonUtils.ExceptionHandler as ExHa
import CommonUtils.progressbar as progBar
from CommonUtils.LoggingManager import Logger

# pyMSsafe modules


class mgftools:
    def __init__(self, hdf5file):
        """
        @brief initialises the mgftools class
        @param hdf5file <string/path>: path for the hdf5 file to analyse
        """
        self.cfgFilters = cfg.parameters['msmsfilters']
        self.maxint = 0
        self.proton = cfg.parameters['general']['proton']
        self.neutron = cfg.parameters['general']['neutron']

        self.hdf5 = hdf5Base(hdf5file)
        self.hdf5.readOpen()

        self.filters = dict(ten=self.tenpercent, zone=self.zonefilter, repion=self.repionfilter,
                            deconv=self.deconvolute, mascot=self.mascot, neutralloss=self.neutrallossfilter,
                            multi=self.multichargefilt, immonium=self.immoniumfilter, none='')
        self.isos = self.hdf5.readTable('/rawdata/isotopes')

        frags = self.getActivationTypes()

        if len(frags) == 1:
            if frags == ['HCD']:
                # hcd olny method
                usefilts = self.cfgFilters['filt_hcd']
                self.remove = self.cfgFilters['rem_hcd']
            else:
                # CID/PQD only method
                usefilts = self.cfgFilters['filt_other']
                self.remove = self.cfgFilters['rem_other']
        else:
            # mixed method
            usefilts = self.cfgFilters['filt_mixed_hcd_other']
            self.remove = self.cfgFilters['rem_mixed_hcd_other']

        if len(self.isos) == 0:
            repionIdx = -1
            for idx in range(len(usefilts)):
                if usefilts[idx] == 'repion':
                    repionIdx = idx
                    break
            if repionIdx != -1:
                usefilts.pop(repionIdx)

        self.usefilts = usefilts
        self.hdf5.close()

    def close(self):
        self.hdf5.close()

    def getActivationTypes(self):

        data = self.hdf5.getDataEqual('/rawdata/parameters', 'parameter', 'Activation Type')
        if len(data) == 0:
            data = self.hdf5.getDataEqual('/rawdata/parameters', 'parameter', 'activation')

        types = {}
        for act in data.flat:
            activation = act['value'].upper()
            types[activation] = types.get(activation, 0) + 1

        return types.keys()

    def export(self, hcdonly=0):
        """
        @brief creates an mgf file for the MS/MS spectra in the hdf5 file
        @param hcdonly <integer>: flag to switch output to only HCD spectra bypassing the normal export filters
        """
        if hcdonly:
            remove = ['CID']
            filters = ['none']
        else:
            remove = self.remove
            filters = self.usefilts

        hdf = self.hdf5
        mgfFile = hdf.filePath.parent.joinpath(hdf.filePath.stem + '.mgf')
        # extra = path(hdf.filepath.splitext()[0] + '.txt')
        # self.fextra = extra.open('w')
        # self.fextra.write('spec_id\tmz\tinten\trel_inten\tion\texp_mz\n')

        mgfOut = open(str(mgfFile), 'w')
        mgfOut.write('#Removed spectra = %s, filtering = %s\n' % (remove, filters))
        spec = 0

        # read parameters from hdf5 file
        try:
            hdf.appendOpen()
            headers = hdf.readTable('/rawdata/msmsheader')
            runTimeEntry = hdf.getDataEqual('/rawdata/parameters', 'parameter', 'MS Run Time (min)')
            if len(runTimeEntry) == 0:
                raise ExHa.MGFprocessingError('MGF Error: Could not find "MS Run Time (min)" parameter in HDF5 file.')
            runtime = runTimeEntry[0]['value']
            units = self.readUnitsOK()

            # add new table for the deconvoluted spectrum data
            hdf.removeTable('/rawdata/deconvions')
            hdf.createTable('rawdata', 'deconvions', 'DeconvIons')
            ident = []
            for frag in units[1]:
                # find all the frag methods to be used in identification
                if 'I' in frag['use']:
                    ident.append(frag['order'])

            logger.log.info('Reading %d spectra from %s' % (len(headers), hdf.filePath.name))
            if 'deconv' in filters:
                deconv = 1
            else:
                deconv = 0

            pBar = progBar.ProgressBar(widgets=progBar.name_widgets, maxval=len(headers), name='Create .mgf').start()
            for idx, h in enumerate(headers):
                if hcdonly:
                    if h['fragmeth'] != 'HCD':
                        continue
                elif not h['order'] in ident:
                    continue
                pBar.update(idx)

                # get spectrum data
                spec = h['spec_id']
                spectrum = hdf.getDataEqual('/rawdata/ions', 'spec_id', spec)
                if deconv:
                    # need extra column for charge information
                    spectrum = self.addChargeColumn(spectrum)

                data = hdf.getDataGeneral('/rawdata/specparams', '(spec_id == %i) & (parameter == "%s")' %
                                          (spec, 'setmass1'))
                setmass = data[0]['value']
                data = hdf.getDataGeneral('/rawdata/specparams', '(spec_id == %i) & (parameter == "%s")' %
                                          (spec, 'frag1'))
                frag = data[0]['value']
                try:
                    self.maxint = max(spectrum['inten'])
                except:
                    self.maxint = 0

                # construct title values list
                rt = '%.3f' % h['rt']
                use = units[1][h['order'] - 1]['use']
                pretitle = ''
                if use == 'IQ':
                    # spec is both ID and Quan so us normal msms ID
                    titles = ['msmsid:F%06d' % h['spec_id']]
                elif use == 'I':
                    if h['quan_spec'] == 0:
                        # no quant data so use spec_id
                        titles = ['msmsid:F%06d' % h['spec_id']]
                    else:
                        # spec is only for ident find the quan spec
                        titles = ['msmsid:F%06d' % h['quan_spec']]
                        pretitle = '#CID=F%06d\n' % h['id_spec']
                elif use == 'Q':
                    titles = ['msmsid:F%06d' % h['quan_spec']]
                    pretitle = '#CID=F%06d\n' % h['id_spec']

                titles.append('rt:' + rt)
                titles.append('survey:S%06d' % h['survey_spec'])
                titles.append('parent:' + setmass)
                titles.append('AnalTime:' + runtime)
                titles.append('Activation:' + frag.upper())

                titleline = 'TITLE=%s\n' % ','.join(titles)

                if h['precmz'] > 0:
                    pepmass = h['precmz']
                elif h['precmz_surv'] > 0:
                    pepmass = h['precmz_surv']
                else:
                    pepmass = h['monomz']

                if pepmass == 0:
                    continue

                for filt in filters:
                    if len(spectrum) > 5 and self.filters[filt]:
                        spectrum = self.filters[filt](h, spectrum)

                # filter for mascot interference
                ionList = []
                if len(spectrum) > 2:
                    mgfOut.write(pretitle)
                    mgfOut.write('BEGIN IONS\n')
                    mgfOut.write(titleline)
                    mgfOut.write('PEPMASS=%f\n' % pepmass)
                    mgfOut.write('CHARGE=%d+\n' % h['charge'])
                    if deconv:
                        for pt in spectrum:
                            if pt['inten'] == 0:
                                continue
                            mgfOut.write('%f  %f  %s\n' % (pt['mz'], pt['inten'], pt['charge']))
                            ionList.append(
                                dict(spec_id=pt['spec_id'], mz=pt['mz'], inten=pt['inten'], charge=pt['charge']))
                    else:
                        for pt in spectrum:
                            if pt['inten'] == 0:
                                continue
                            mgfOut.write('%f  %f\n' % (pt['mz'], pt['inten']))
                            ionList.append(dict(spec_id=pt['spec_id'], mz=pt['mz'], inten=pt['inten']))
                    mgfOut.write('END IONS\n\n')
                if len(ionList) > 0:
                    hdf.appendRows('/rawdata/deconvions', ionList)

            pBar.finish()

        except ExHa.MGFprocessingError, czEx:
            if spec:
                ExHa.addContext(czEx, 'Raised whist processing spectrum %i' % spec)
            raise
        except Exception as genEx:
            ExHa.reformatException(genEx)
            if spec:
                ExHa.addContext(genEx, 'Raised whist processing spectrum %i' % spec)
            raise

        mgfOut.close()
        hdf.indexTable('/rawdata/deconvions', ['spec_id'])
        hdf.close()
        return {'code': 0}
        # self.fextra.close()

    def readUnitsOK(self):
        data = self.hdf5.getDataGeneral('/rawdata/units', '(unit > 0)')
        units = {}
        for d in data:
            units.setdefault(d['unit'], [])
            unit = dict(order=d['order'], activation=d['activation'], acttime=d['acttime'],
                        energy=d['energy'], isolation=d['isolation'], lowmz=d['lowmz'],
                        resolution=d['resolution'], use=d['use'], scanevents=eval(d['scanevents']),
                        colenergysteps=d['colenergysteps'], colenergywidth=d['colenergywidth'])
            units[d['unit']].append(unit.copy())

        return units

    @staticmethod
    def addChargeColumn(spectrum):
        charge = np.ndarray(len(spectrum), dtype=[('spec_id', int), ('mz', float), ('inten', float), ('charge', 'S9')])
        for i, ion in enumerate(spectrum):
            charge.put(i, (ion['spec_id'], ion['mz'], ion['inten'], '0'))

        return charge

    # noinspection PyUnusedLocal
    @staticmethod
    def tenpercent(header, spectrum):
        tenperc = len(spectrum) / 10
        spectrum.sort(order='inten')
        tmp = spectrum[tenperc:]
        tmp.sort(order='mz')

        return tmp

    def zonefilter(self, header, ions):
        """
        @brief does the Matheus Mann filtering top x ions in a chunk of the spectrum (50 / 100 m/z)
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        """
        # first remove 1+ isotope data
        ions = self.doDeisotope(ions)
        numions = len(ions)
        if numions < 10:
            # print 'only %d ions present - no filtering done' % numions
            return ions

        bindic = {}
        if header['charge'] <= 2:
            chunk = 100
        else:
            chunk = 50
        ions.sort(order='inten')
        for i in range(numions - 1, -1, -1):
            if ions[i]['inten'] == 0:
                break

            mzBin = round(ions[i]['mz'] / chunk, 0) * chunk
            num = len(bindic.get(mzBin, ''))
            if num <= 4:
                bindic.setdefault(mzBin, []).append(i)
            else:
                # not to be included so set inten to zero
                ions[i]['inten'] = 0

        # remove all ions with zero intensity
        ions.sort(order='inten')
        i = 0
        for i in range(numions):
            if ions[i]['inten'] > 0:
                break

        ions = ions[i:]

        # set back to mass ordered array
        ions.sort(order='mz')
        maxmass = (header['monomz'] - self.proton) * header['charge'] - 50

        last = len(ions) - 1
        while ions[last]['mz'] > maxmass:
            last -= 1

        return ions[:last + 1]

    def neutrallossfilter(self, header, ions):
        """
        @brief filters the neutral loss ions  from the deconvoluted precursor ion
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        """
        # proton = self.cfg.proton
        neutrals = self.cfgFilters['neutrals']
        if self.usefilts.index('deconv') < self.usefilts.index('neutralloss'):
            # data has been deconvoluted before the NL filter is applied
            # need to use the deconvoluted precursor mass
            prec = (header['precmz'] - self.proton) * header['charge'] + self.proton
        else:
            prec = header['precmz']
            for idx, p in enumerate(neutrals):
                neutrals[idx] = p / header['charge']

        nlions = []
        for delta in neutrals:
            nlions.append(prec - delta)

        filtered, removed = self.largeionfilter(nlions, ions, self.cfgFilters['tolmda'], self.cfgFilters['tolppm'])

        return filtered

    # noinspection PyUnusedLocal
    def immoniumfilter(self, header, ions):
        """
        @brief removes reporter ions from the MS/MS spectrum
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        """
        masses = self.cfgFilters['immonium']
        filtered, removed = self.smallionfilter(masses, ions, self.cfgFilters['tolmda'], self.cfgFilters['tolppm'])

        return filtered

    # noinspection PyUnusedLocal
    def repionfilter(self, header, ions):
        """
        @brief removes reporter ions from the MS/MS spectrum
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        """
        masses = []
        for i in self.isos:
            masses.append(i['mz'])

        filtered, removed = self.smallionfilter(masses, ions, self.cfgFilters['tolmda'], self.cfgFilters['tolppm'])

        return filtered

    def deconvolute(self, header, ions):
        """
        @brief deconvolutes the MS/MS spectrum
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        """
        proton = self.proton

        # remove 13C ions and identify charge states
        deiso, removed = self.deisotope(header, ions)
        maxcharge = 1
        if removed:
            # only do the deconvolution if data was removed ie isotopes were detected
            for ion in deiso:
                cge = float(ion['charge'])
                if cge > 1:
                    deconv = (ion['mz'] - proton) * cge + proton
                    ion['mz'] = deconv
                    if cge > maxcharge:
                        maxcharge = cge

        deiso.sort(order='mz')
        # 1.5 is abitrary value to ensure the missincorp of label is never seen as an ion
        if header['charge'] == 0:
            maxmass = (header['monomz'] - proton) * maxcharge - 1.50
        else:
            maxmass = (header['monomz'] - proton) * header['charge'] - 1.50

        last = len(deiso) - 1
        while deiso[last]['mz'] > maxmass and last >= 0:
            last -= 1

        return deiso[:last + 1]

    def deisotope(self, header, ions):
        """
        @brief removes 13C isotopes from the MS/MS spectrum
        @param header <np.array>: containing the spectrum header data
        @param ions <np.array>: containing the ion data
        @return filtered <np.array>: filtered data
        @return removed <integer>: the number of ions removed from the spectrum
        """
        neutron = self.neutron
        removed = 0

        # calculate the possible deltas for the precursor charge
        neutrondelta = [(neutron / c, c) for c in range(header['charge'], 0, -1)]

        # loop backwards through the data
        for c13 in range(len(ions) - 1, 0, -1):
            mz13c = ions[c13]['mz']
            if mz13c < 132:
                break
            c12 = c13
            delta = 0
            # now scan back for a prior ion
            while delta < 2 and c12 > 0:
                c12 -= 1
                delta = mz13c - ions[c12]['mz']
                for nd in neutrondelta:
                    if abs(delta - nd[0]) < 0.02 and ions[c12]['inten'] > ions[c13]['inten'] / 2:
                        # ion matches in mass and intensity
                        cge = str(nd[1])
                        ions[c12]['charge'] = cge
                        ions[c13]['charge'] = cge
                        ions[c13]['inten'] = 0
                        removed += 1

        return ions, removed

    @staticmethod
    def doDeisotope(ions):
        """
        @brief simple method to remove 13C isotopes
        @param ions <np.array>: containing the ion data
        @return ions <np.array>: data with 13C peak intensities set to zero
        """
        for c13 in range(len(ions) - 1, -1, -1):
            if 133 > ions[c13]['mz'] > 113:
                ions[c13]['inten'] = 0
            for c12 in range(c13 - 1, -1, -1):
                delta = ions[c13]['mz'] - ions[c12]['mz']
                if delta > 2:
                    break
                elif abs(delta - 1.0033548) < 0.3 and ions[c12]['inten'] > ions[c13]['inten']:
                    ions[c13]['inten'] = 0

        return ions

    # noinspection PyUnusedLocal
    @staticmethod
    def mascot(header, ions):
        """
        @brief method to check data for ions that would be reassigned by mascot
        @param ions <np.array>: containing the ion data
        @return ions <np.array>: data with all ions independent of mascot reassignment
        """
        length = len(ions)

        for i, ion in enumerate(ions):
            j = i + 1
            if ion['inten'] == 0 or j >= length:
                break  # don't process ions with zero intensity
            delta = ions[j]['mz'] - ion['mz']
            tmp = [(ion['mz'], ion['inten'], ion['charge'])]
            while delta < 0.01:
                if ions[j]['inten'] > 0:
                    tmp.append((ions[j]['mz'], ions[j]['inten'], ions[j]['charge']))
                    ions[j]['inten'] = 0
                j += 1
                if j >= length:
                    break
                delta = ions[j]['mz'] - ion['mz']

            if len(tmp) > 1:
                # found ions to filter so need to recalculate the m/z and inten
                sumprod = 0.0
                suminten = 0.0
                for t in tmp:
                    sumprod += t[0] * t[1]
                    suminten += t[1]
                newmz = sumprod / suminten

                ions.put(i, (ion['spec_id'], newmz, suminten, '0'))

        return ions

    # noinspection PyUnusedLocal
    def multichargefilt(self, header, ions):
        """
        @brief method to check data for ions that would be reassigned by mascot
        @param ions <np.array>: containing the ion data
        @return ions <np.array>: data with all ions independent of mascot reassignment
        """
        length = len(ions)
        ppm = self.cfgFilters['tolppm']
        mda = self.cfgFilters['tolmda']

        for i, ion in enumerate(ions):
            j = i + 1
            # only process ions with valid m/z and intensities and charges
            if ion['inten'] == 0 or ion['mz'] == 0:
                continue  # or ion['charge'] == '0.0'
            if j >= length:
                break  # stop when the end of the array is reached
            delta = ions[j]['mz'] - ion['mz']
            error = max(mda, ion['mz'] * ppm)
            tmp = [(ion['mz'], ion['inten'], ion['charge'])]
            while delta <= error:
                if ions[j]['charge'] != '0.0' or ion['charge'] != '0.0':
                    tmp.append((ions[j]['mz'], ions[j]['inten'], ions[j]['charge']))
                    ions[j]['mz'] = 0
                j += 1
                if j >= length:
                    break
                delta = ions[j]['mz'] - tmp[-1][0]

            if len(tmp) > 1:
                sumprod = 0.0
                suminten = 0.0
                cge = []
                for t in tmp:
                    sumprod += t[0] * t[1]
                    suminten += t[1]
                    cge.append(t[2])
                newmz = sumprod / suminten
                ions.put(i, (ion['spec_id'], newmz, suminten, ','.join(cge)))
        ions.sort(order='mz')
        i = 0
        while ions[i]['mz'] == 0:
            i += 1

        return ions[i:]

    @staticmethod
    def largeionfilter(masses, spec, mda, ppm):
        """
        @brief general ion filter, works from the high mass end of the ions array
        @param masses <list>: masses to be removed from the spectrum
        @param spec <np.array>: containing the ion data
        @return filtered <np.array>, removed <list>: filtered data and data from the ions removed from the spectrum
        """

        d = 0
        error = max(mda, masses[d] * ppm)

        i = len(spec) - 1
        removed = []
        while i >= 0:
            ion = spec[i]
            delta = ion['mz'] - masses[d]
            if abs(delta) <= error:
                # ion needs to be filtered
                removed.append((ion['spec_id'], ion['mz'], ion['inten'], d))
                ion['mz'] = 0
                i -= 1
            elif delta < - error:
                # below range of this ion
                d += 1
                if d >= len(masses):
                    break
                error = max(mda, masses[d] * ppm)
            else:
                i -= 1

        # ions were removed so filter them from the dataset
        if removed:
            spec.sort(order='mz')
            spec = spec[len(removed):]

        return spec, removed

    @staticmethod
    def smallionfilter(masses, spec, mda, ppm):
        """
        @brief general ion filter, works from the low mass end of the ions array
        @param masses <list>: masses to be removed from the spectrum
        @param spec <np.array>: containing the ion data
        @return filtered <np.array>, removed <list>: filtered data and data from the ions removed from the spectrum
        """

        d = 0
        error = max(mda, masses[d] * ppm)

        i = 0
        removed = []
        lenspec = len(spec)
        while i < lenspec:
            ion = spec[i]
            delta = ion['mz'] - masses[d]
            if abs(delta) <= error:
                # ion needs to be filtered
                removed.append((ion['spec_id'], ion['mz'], ion['inten'], d))
                ion['mz'] = 0
                i += 1
            elif delta > error:
                # below range of this ion
                d += 1
                if d >= len(masses):
                    break
                error = max(mda, masses[d] * ppm)
            else:
                i += 1

        # ions were removed so filter them from the dataset
        if removed:
            spec.sort(order='mz')
            spec = spec[len(removed):]

        return spec, removed


# ########### MAIN ###############

if __name__ == '__main__':
    f = 0
    try:
        configPath = './mgf.cfg'
        cfg = pyMSsafeConfigManager(configPath)
        ret = cfg.evaluateCommandLineArgs(sys.argv)
        cfg.scalePpmMda()

        dataDir = cfg.parameters['runtime']['datadir']
        fileFilter = cfg.parameters['runtime']['filefilter']
        hcdOnly = cfg.parameters['general']['hcdonly']

        logParam = cfg.parameters['logging']
        logPath = Path(dataDir.joinpath(logParam['logdir']))
        if not logPath.exists():
            logPath.mkdir(parents=True)
        logFile = logPath.joinpath(logParam['logfile'])

        logger = Logger(logFile, logParam['loglevel'], logParam['screenlevel'], False)
        logger.setLog('mgfCreation')

        logger.log.info('Analysing: %s' % str(dataDir.joinpath(fileFilter)))
        if hcdOnly:
            logger.log.info('Exporting HCD data only')

        # for f in dataDir.files(fileFilter):
        for f in dataDir.glob(fileFilter):
            if not f.is_file():
                # skip any directories
                continue

            # if f.name[:4] in ['6528', '1814', '2032']: continue
            mgf = mgftools(f)
            logger.log.info('Filename:     %s' % f.name)
            if hcdOnly:
                logger.log.info('Export HCD data only')
            else:
                logger.log.info('Using filters: %s' % str(mgf.usefilts))
            rtn = mgf.export(hcdOnly)
            mgf.close()
        if f == 0:
            raise ExHa.MGFprocessingError('no files found for: %s' % str(dataDir / fileFilter))

    except ExHa.UsageError as useEx:
        ExHa.reformatException(useEx)
        logger.log.info(useEx.context)
    except Exception, genEx:
        ExHa.reformatException(genEx)
        if f:
            ExHa.exportError2File(genEx, f.parent.joinpath(f.stem + '.error'))
        else:
            ExHa.exportError2File(genEx, dataDir.joinpath('errors.error'))
        logger.log.info(ExHa.multiLineRepr(genEx))

    logger.log.info('finished')
