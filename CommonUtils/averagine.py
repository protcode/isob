from pathlib import Path
from ChemicalDefinitions import ChemicalDefinitions
from IsotopeCalculator import IsotopeCalculator


class AveragineModel:
    """
    @brief class to manage theoretical isotope distributions of peptides
    """

    def __init__(self, dataDir, maxRepeats, minInten, quanMeth, numLabels, cfg, overwrite=0):
        """
        @brief initiates the averagine isotopic abundance model
        @param isofile <path object>: linking to the idotope data file
        @param maxrepeats <integer>: maximum number of averagines to calculate the isotopic abundances for
        @param mininten <float>: minimum fraction of the normalised intensity for an isotope to be used
        @param quanMeth <string>: the quantification method used to select the correct isofile
        """
        self.maxRepeats = maxRepeats
        self.minInten = minInten
        self.c12 = None
        self.isolist = None
        self.numlabels = numLabels
        self.CD = ChemicalDefinitions()
        self.IC = IsotopeCalculator(cfg, self.CD)
        self.cfg = cfg

        self.calcIsotopeDistribution(quanMeth)

        isof = dataDir.joinpath('averagine_%s.txt' % quanMeth)
        #
        # if not isof.exists() or overwrite:
        #     self.calcIsotopeDistribution(quanMeth)
        #     self.savedata(isof, quanMeth)
        # else:
        #     self.loadfile(isof)

    def loadfile(self, isof):
        """
        @brief loads the isotope ratios from file
        @param isof <object>: path object for the isoratio file
        """
        fin = open(str(isof), 'r')
        line = fin.readline()
        splt = line[:-1].split()
        limits = (float(splt[0]), int(splt[1]))
        self.c12 = int(splt[2])
        self.maxRepeats = limits[1]
        isos = []
        for j in range(limits[1]):
            line = fin.readline()
            data = line[:-1].split('\t')
            numreps = int(data[0])
            mass = float(data[1])
            comp = data[2]
            ratios = [float(x) for x in data[3:]]
            ratiosDict = {}
            for idx, inten in enumerate(ratios):
                ratiosDict[idx - self.c12] = inten
            isos.append(dict(monomass=mass, numreps=numreps, comp=comp, inten=ratios[:], RelInts=ratiosDict.copy(),
                             sumTheo=sum(ratios)))
        self.isolist = isos
        fin.close()

    def findc13inten(self, mass):
        """
        @brief finds the closest averagine model to the mass and returns the intensity of the +1 13C isotope
        @param mass <float>: the mass to find the averagine for
        @return inten <float>: the averagine 13C intensity
        """

        isos = self.findIsotopeSet(mass)
        return isos['inten'][self.c12 + 1]

    def findIsotopeSet(self, mass, remedge=0):
        """
        @brief finds the closest averagine model to the mass
        @param mass <float>: the mass to find the averagine for
        @param remedge <integer>: flag to remove isotopes below the 12C, 1 removes below 12C 2 also removes high isos
                according to mass
        @return iso <tuple>: the averagine isotope model (m/z <float>, num isos <integer>, isotope ratios
                <list of floats>)
        """
        a = 0
        while a < self.maxRepeats:
            if mass < self.isolist[a]['monomass']:
                break
            a += 1
        if a == self.maxRepeats:
            iso = self.isolist[-1]
        elif abs(mass - self.isolist[a]['monomass']) <= abs(mass - self.isolist[a - 1]['monomass']):
            iso = self.isolist[a]
        else:
            iso = self.isolist[a - 1]

        # calculate the isotpe intensities if needed
        if 'inten' not in iso:
            # need to calculate the isotope data
            isoDist = self.IC.filterRelIntensities(self.IC.calcIsotopeDistributionFromElemComp(iso['elements']))
            iso.update(isoDist)

            # find the isotope with the greatest intensity
            maxpos = 0
            maxint = 0
            for i in range(len(iso['inten'])):
                if iso['inten'][i] > maxint:
                    maxint = iso['inten'][i]
                    maxpos = i
            iso['max'] = maxpos - iso['c12inten']

            if remedge >= 1:
                iso['inten'] = iso['inten'][iso['c12inten']:]
                maxpos -= iso['c12inten']

            if remedge == 2:
                edge = maxpos + self.cfg.parameters['deisotoping']['num_past_max'] + 1
                if edge > 4:
                    iso['inten'] = iso['inten'][:edge]
                else:
                    iso['inten'] = iso['inten'][:4]

        return iso

    def calcIsotopeDistribution(self, quanMeth, maxRepeats=0, minInten=0, overwrite=0):
        """
        @brief calculates the averagine isotopic abundances
        @param quanMeth <string>: the quantification method used to select the correct isofile
        """

        elem = self.CD.averagine_composition

        if maxRepeats:
            self.maxRepeats = maxRepeats
        if minInten:
            self.minInten = minInten

        self.isolist = []
        # iso = dict(monomass=0.0, avgmass=0.0, inten=[], norm=[], c12index=0, C=0, H=0, N=0, O=0, S=0, C13=0, N15=0)

        # start with the labels if any
        elemComp = self.QuantificationLables(quanMeth)

        for rep in range(1, self.maxRepeats + 1):
            # bild up the isotopic composition
            for e in elem:
                try:
                    elemComp[e] += elem[e]
                except KeyError:
                    elemComp[e] = elem[e]

            # create a composition for this repeat
            repComp = {}
            for e in elemComp:
                repComp[e] = int(elemComp[e] + 0.5)

            averagineData = self.IC.calcMassesFromElemComp(repComp)
            averagineData.update(self.IC.calcIsotopeDistributionFromElemComp(repComp))
            averagineData = self.IC.filterRelIntensities(averagineData)
            averagineData['comp'] = self.elements2string(repComp)
            averagineData['numreps'] = rep
            averagineData['RelInts'] = {}
            averagineData['sumTheo'] = sum(averagineData['inten'])
            for idx, inten in enumerate(averagineData['inten']):
                averagineData['RelInts'][idx - averagineData['c12index']] = inten
            self.isolist.append(averagineData.copy())

        return 1

    def savedata(self, isof, quanMeth):
        """
        @brief saves the calculated isotope data to file
        @param isof <object>: path object for the isoratio file
        @param quanMeth <string>: the quantification method used to select the correct isofile
        """
        fout = open(str(isof), 'w')
        c12 = self.isolist[0]['c12index']
        fout.write('%.3f\t%d\t%d\t%s\t%d\n' % (self.minInten, self.maxRepeats, c12, quanMeth, self.numlabels))

        for rep in self.isolist:
            # ouput the normal element data
            line = '%i\t%.6f\t%s' % (rep['numreps'], rep['monomass'], self.elements2string(rep['elements']))

            # output the labeled elemental data
            for label in ['2H', '13C', '15N', '18O']:
                if label in rep['elements'] and rep['elements'][label] > 0:
                    line += ' %s%d' % (label, rep['elements'][label])

            # output the intensity data
            if rep['c12inten'] < c12:
                # fewer isotopes before the c12 then add
                line += '\t0.0' * (c12 - rep['c12inten'])

            for i in rep['inten']:
                line += '\t%.4f' % i

            fout.write(line + '\n')

        fout.close()

    def elements2string(self, elements):
        elemStr = 'C%d H%d N%d O%d S%d' % (elements['C'], elements['H'], elements['N'], elements['O'], elements['S'])

        for label in ['2H', '13C', '15N', '18O']:
            if label in elements and elements[label] > 0:
                elemStr += ' %s%d' % (label, elements[label])

        return elemStr

    def QuantificationLables(self, quanMeth):
        if quanMeth in self.CD.heavy_isotopes:
            labelComposition = self.CD.heavy_isotopes[quanMeth].copy()
        else:
            labelComposition = self.CD.heavy_isotopes['none'].copy()

        for e in labelComposition:
            labelComposition[e] = int(labelComposition[e] * self.numlabels + 0.5)

        return labelComposition

# ########### MAIN ###############
if __name__ == '__main__':

    class config:
        def __init__(self):
            self.parameters = dict(isotopecalc=dict(calcmininten=0.0001, reportmininten=0.01, minmissincorp=0))
    cfg = config()

    def findData():
        for mz in [1000, 5000, 20000]:
            isos = am.findIsotopeSet(mz)
            c13 = am.findc13inten(mz)
            print '\t%i\t%i\t%.4f\t%s\t%.4f\t%s' % (mz, isos['numreps'], isos['monomass'], isos['comp'], c13,
                                                    isos['inten'][:4])

    for quanMeth in ['none', 'tmt', 'itraq', '8traq']:
        print 'Processing quantification method =', quanMeth
        print 'creating new file'
        am = AveragineModel(Path('.'), 200, 0.01, quanMeth, 2, cfg, overwrite=1)
        findData()
        print 'reading from file'
        am = AveragineModel(Path('.'), 200, 0.01, quanMeth, 2, cfg, overwrite=0)
        findData()
