#!/usr/bin/env python

'''
basic script to copy files
'''

import sys
import glob
import os
import shutil
from pathlib import Path

src = Path(sys.argv[1])  # source FOLDER
dst = Path(sys.argv[2])  # source file called dst/contig.fasta

try:
    segs = str(sys.argv[3]).lower()
except IndexError:
    sys.exit("Did you specify argument 3 - subset or all?")

print('\nConcatenating files from: {source} to {destination} for {segments} segments \n'.format(
    source=src, destination=dst, segments=segs))

if segs in ['all', 'subset', 'rsv']:
    if segs == 'all':
        seg_to_keep = ['HA', 'NA', 'MP', 'NS', 'NP', 'PA', 'PB1', 'PB2']
    elif segs == 'subset':
        seg_to_keep = ['HA', 'NA', 'MP']
    elif segs == 'rsv':
        seg_to_keep = ['rsv']
else:
    sys.exit("Error - subset (argument #3) must be all or subset")


filenames = glob.iglob(os.path.join(src, "*.fasta"))
with open(dst, 'w') as outfile:
    for fname in filenames:
        if Path(fname).stem in seg_to_keep:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            outfile.write('\n')


exit(0)
