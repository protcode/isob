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

This file is the interface to the quantification methods defined in
QuantMethod.cfg, providing access and lookup functions for the methods.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob
"""

# python libraries
import copy
import sys
from pathlib import Path

# CommonUtils libraries
sys.path.insert(0, '..')
from CommonUtils.ConfigManager import ConfigManager

# application libraries


class QuantMethods:
    def __init__(self):
        modulePath = Path(__file__).parent
        for pth in modulePath.parents:
            configPath = pth.joinpath('QuantMethod.cfg')
            if configPath.exists():
                break

        cfg = ConfigManager(str(configPath))
        self.methodsByName = {}
        self.methodsByID = {}
        for s in cfg.parameters:
            meth = cfg.parameters[s]
            method = {'quantmasses': {}}
            for key in meth:
                if key in ['__name__']:
                    continue
                elif key in ['method_id', 'batch_id', 'num_isotopes']:
                    method[key] = int(meth[key])
                elif key in ['hyperplexed', 'heavy_isotopes']:
                    method[key] = eval(meth[key])
                elif key in ['xic_only']:
                    if meth[key].lower() == 'yes' or meth[key] == '1' or meth[key] == 'True':
                        method[key] = True
                    else:
                        method[key] = False
                else:
                    try:
                        iso = int(key)
                        method['quantmasses'][iso] = eval(meth[key])
                    except ValueError:
                        method[key] = meth[key]

            self.methodsByID[method['method_id']] = copy.deepcopy(method)
            self.methodsByName[method['method_name'].upper()] = copy.deepcopy(method)

    def getMethodByID(self, method_id):
        if method_id in self.methodsByID:
            return self.methodsByID[method_id]
        else:
            return None

    def getMethodByName(self, method_name):
        if method_name in self.methodsByName:
            return self.methodsByName[method_name]
        else:
            return None

    def getMethodByIsotope(self, isotope):
        for id in self.methodsByID:
            meth = self.methodsByID[id]
            if isotope in meth['quantmasses']:
                return meth
        return None

    def extractCorrectionFactors(self, method):
        correction = {}
        for id in method['quantmasses']:
            if 'corrections' in method['quantmasses'][id][0]:
                correction[id] = method['quantmasses'][id][0]['corrections']
            else:
                correction[id] = {}
        return correction


# ########### MAIN ###############
if __name__ == '__main__':
    qmh = QuantMethods()

    print qmh.getMethodByName('TMT10')
    print
    print qmh.getMethodByID(12)
    print
    print qmh.getMethodByIsotope(78)
    print
    print qmh.extractCorrectionFactors(qmh.getMethodByID(12))
