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

This file controls the parsing of the unimod XML data embedded in
Mascot .dat files.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/
"""

from xml.sax import make_parser, ContentHandler


class UnimodXMLParser(ContentHandler):
    def __init__(self):
        self.inAA = False
        self.inMod = False
        self.inDelta = False
        self.inSpecificity = False
        self.inNeutLoss = False
        pass

    def startElement(self, name, attrs):
        """
        Deals with the finding a new xml tag.
        """

        if name == 'umod:unimod':
            self.unimod = attrs._attrs.copy()
        elif name == 'umod:elements':
            self.elements = {}
            pass
        elif name == 'umod:elem':
            self.elements[attrs['title']] = attrs._attrs.copy()
            pass
        elif name == 'umod:modifications':
            self.mods = {}
            pass
        elif name == 'umod:mod':
            self.inMod = True
            self.mod = dict(title=attrs['title'], full_name=attrs['full_name'], specificity={}, da_delta={})
            pass
        elif name == 'umod:amino_acids':
            self.aminoacids = {}
            pass
        elif name == 'umod:aa':
            self.inAA = True
            self.amino = attrs._attrs.copy()
            self.amino['elements'] = {}
        elif name == 'umod:element':
            # element used in several locations
            if self.inAA:
                addto = self.amino['elements']
            elif self.inDelta:
                addto = self.delta['elements']
            elif self.inSpecificity:
                addto = self.nl['elements']

            addto[attrs['symbol']] = int(attrs['number'])
        elif name == 'umod:delta':
            self.inDelta = True
            self.delta = attrs._attrs.copy()
            self.delta['elements'] = {}
        elif name == 'umod:specificity':
            self.inSpecificity = True
            self.specif = attrs._attrs.copy()
            self.specif['neutral_loss'] = []
        elif name == 'umod:NeutralLoss':
            self.inNeutLoss = True
            self.nl = attrs._attrs.copy()
            self.nl['elements'] = {}

    def endElement(self, name):
        if name == 'umod:aa':
            aa = self.amino
            if aa['title'] != '-':
                self.aminoacids[aa['title']] = aa.copy()
            self.inAA = False
        elif name == 'umod:mod':
            self.inMod = False
            self.mods[self.mod['title']] = self.mod.copy()
        elif name == 'umod:delta':
            self.inDelta = False
            self.mod['delta'] = self.delta.copy()
        elif name == 'umod:specificity':
            self.inSpecificity = False
            self.mod['specificity'][self.specif['spec_group']] = self.specif.copy()
        elif name == 'umod:NeutralLoss':
            self.inNeutLoss = False
            if not self.nl['composition'] in ['0', '']:
                self.specif['neutral_loss'].append(self.nl.copy())

if __name__ == "__main__":
    from pathlib import Path

    file = Path('./data/tmpold.xml')
    tmpxml = file.open('r')
    xml = UnimodXMLParser()
    saxparser = make_parser()
    saxparser.setContentHandler(xml)
    saxparser.parse(tmpxml)
    tmpxml.close()
