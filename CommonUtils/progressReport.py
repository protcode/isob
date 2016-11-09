__author__ = 'gs625046'


class progressReport:

    def __init__(self, maxCount, filename, process, unit):
        self.maxCount = float(maxCount)
        self.reportStep = 0.2
        self.nextReport = self.reportStep
        self.counter = 0
        self.baseProgressString = '    %-20s from %-30s: ' % (process, filename)
        print self.baseProgressString + 'Starting %i %s' % (maxCount, unit)

    def next(self):
        self.counter += 1
        self.report(self.counter)

    def report(self, counter):
        progress = counter / self.maxCount
        if progress >= self.nextReport:
            print self.baseProgressString + '%i%% completed' % int(self.nextReport * 100)
            self.nextReport += self.reportStep

    def endReport(self):
        print self.baseProgressString + 'Finished'
