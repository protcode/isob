###############################################################
#This config file is part of the isobarQuant package,
#written by Toby Mathieson and Gavain Sweetman
#(c) 2016 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
#69117, Heidelberg, Germany.

#The isobarQuant package processes data from
#.raw files acquired on Thermo Scientific Orbitrap / QExactive / Fusion
#instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
#It creates an .hdf5 file into which are later parsed the results from
#Mascot searches. From these files protein groups are inferred and quantified.

#This config file is used by proteinInference.py

#isobarQuant is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#A copy of the license should have been part of the
#download. Alternatively it can be obtained here:
#https://github.com/protcode/isob
###############################################################

# this should be called from mascot daemon, external processes tab, section 'After each search'
# python c:\isobarQuant\CommonUtils\getDatFiles.py <resultfilepath> <resultfilename>
#           http://mymascotserver/mascot/cgi/ c:/isobarQData/vechicle_1
# resultfilepath ==> ../data/20160101/F123456.dat
# resultfilename ==> F123456.dat

import sys

from pathlib import Path
import urllib

resultfilepath = sys.argv[1]
resultfilename = sys.argv[2]
mascotbaseURL = sys.argv[3].strip()
localdirectory = sys.argv[4].strip()


mascotURL = mascotbaseURL + '/export_dat_2.pl?do_export=1&export_format=MascotDAT&file=%s'
savelocation = str(Path(localdirectory) / Path(resultfilename))


testfile = urllib.URLopener()
testfile.retrieve(mascotURL % resultfilepath, savelocation)
