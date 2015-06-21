"""This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This small script simply extracts the gene names and descriptions from
Uniprot .fasta files, writing them to a file called essentialuniprotdata.txt in
the current directory.

The first command line argument is the location of the .fasta file.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here:
https://github.com/protcode/isob/
"""

import sys
import re

rg_gene = re.compile('GN=(\w+)')

fastafilelocation = sys.argv[1]

fh = open(fastafilelocation, 'r')

id2data = {}
for line in fh:
    if line.startswith('>'):
        items = line.strip().split('|')
        proteinid = items[1]

        description = items[2]
        if rg_gene.search(description):
            gene_name = rg_gene.search(description).group(1)
            if proteinid.startswith('DD') or proteinid.startswith('###REV###') or proteinid.startswith('###RND###'):
                gene_name = '###%s###' % gene_name
        else:
            gene_name = 'NA'
        id2data[proteinid] = (description, gene_name)
ofh = open('essentialuniprotdata.txt', 'w')
for proteinid, data in id2data.iteritems():

    ofh.write('%s\t%s\t%s\n' % (proteinid, data[0], data[1]))

ofh.close()
