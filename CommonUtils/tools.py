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

This file provides an number of miscellaneous tools.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

import time


class Jobcontainer:
    """
    @brief hold information needed for running jobs.
    """

    def __init__(self, source, jobid):
        # create an empty job
        self.srcpth = source
        self.dstpth = source.parent
        self.config = ''
        self.job_id = jobid
        self.args = dict(maxspec=0)
        self.mode = ''


class Stopwatch(object):
    """
    @brief class that deals with timing events in the code operation
    """

    def __init__(self):
        self.start = time.clock()
        self.events = []

    def rec(self, name):
        """
        @brief records a time point with named identifier
        @param name <string>: identifier for the timed event
        """
        self.events.append((name, time.clock()))

    def stop(self):
        """
        @brief ends the timing events
        """
        self.events.append(('runtime', time.clock()))
        self.result = []
        for i in range(len(self.events)):
            x = self.events[i]
            if i == 0:
                t_last = self.start
            else:
                t_last = self.events[i - 1][1]
            self.result.append(dict(name=x[0], split='%.3f' % (x[1] - self.start), lap='%.3f' % (x[1] - t_last),
                                    lapf=x[1] - t_last))
        return self.result

    def oneLineFormat(self):
        out = []
        for i in self.result:
            out.append('%s = %s s' % (i['name'], i['lap']))

        if i:
            out[-1] = '%s = %s s' % (i['name'], i['split'])
        return ', '.join(out)

    def format(self):
        """
        @brief formats the event for screen display
        """
        out = []
        res = 0
        nameSpace = 0
        timeSpace = 0

        for res in self.result:
            if len(res['name']) > nameSpace:
                nameSpace = len(res['name'])

        timeSpace = len(self.result[-1]['split'])

        strFormat = '%-' + str(nameSpace) + 's = %' + str(timeSpace) + 's s'

        for res in self.result:
            out.append(strFormat % (res['name'], res['lap']))

        if res:
            out[-1] = strFormat % (res['name'], res['split'])
        return out

    def write(self):
        """
        @brief writes the events to screen
        """
        t = '\n'.join(self.format())
        print '\n'.join(self.format())


class IO(object):

    @staticmethod
    def numpyToCSV(inputArray, outPath, sep='\t'):
        fout = open(outPath, 'w')

        names = sep.join(inputArray.dtype.names)
        fout.write(names + '\n')

        for row in inputArray:
            rowData = [str(x) for x in row]
            data = sep.join(rowData)
            fout.write(data + '\n')

        fout.close()


class Misc(object):

    @staticmethod
    def uniqifyList(lst):
        """
        @brief removes duplicates from a list while preserving the order of elements in the input list
        @param lst <list>: input list
        @return noDupes <list>: list without duplicates in the same order as the input list
        """
        noDupes = []
        [noDupes.append(i) for i in lst if not noDupes.count(i)]
        return noDupes
