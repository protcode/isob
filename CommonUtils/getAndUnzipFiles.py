__author__ = 'tm224714'

import sys
import shutil
import subprocess
import os
from pathlib import Path

startpnamedirectory = Path(sys.argv[1])
destinationdirectory = Path(sys.argv[2])
if not destinationdirectory.exists():
    destinationdirectory.mkdir()
p = [x for x in startpnamedirectory.glob('*.rar')]
filelist = list()

for file in p:
    if file.name.find('BLANK') > -1 or file.name.find('BSA') > -1:
        continue
    else:
        print 'file', file.name
        copiedfile = destinationdirectory / Path(file.name)
        shutil.copy(str(file), str(copiedfile))
        filelist.append(copiedfile)

for file in filelist:
    cmdline = 'UnRAR.exe -y e %s %s ' % (str(file), str(file.parent))
    print cmdline
    os.system(cmdline)

for file in filelist:
    # and delete .rar file.
    file.unlink()
