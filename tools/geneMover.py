#!/usr/bin/env python

'''
basic script to copy files
'''

import sys
import glob
import os
import shutil
from pathlib import Path

if len(sys.argv) == 1:
    sys.exit('''
    geneMover.py - Combines fasta files from a given directory

    Usage: 
        geneMover.py sourceFolder destination_file [subset, all, rsv]
        geneMover.py output/assembly/ output/assembly/rename/ rsv
    ''')

src = Path(sys.argv[1])  # source folder
dst = Path(sys.argv[2])  # source file called dst/contig.fasta

try:
    segs = str(sys.argv[3]).lower()
except IndexError:
    sys.exit("Did you specify argument 3 - subset, all, rsv?")

print('\nConcatenating files from: {source} to {destination} for {segments} segments \n'.format(
    source=src, destination=dst, segments=segs))

if segs in ['all', 'subset', 'rsv']:
    if segs == 'all':
        seg_to_keep = ['HA', 'NA', 'MP', 'NS', 'NP', 'PA', 'PB1', 'PB2']
    elif segs == 'subset':
        seg_to_keep = ['HA', 'NA', 'MP']
    elif segs == 'rsv':
        seg_to_keep = ['RSV']
else:
    sys.exit("Error - subset (argument #3) must be: all, subset, rsv")


filenames = glob.iglob(os.path.join(src, "*.fasta"))
with open(dst, 'w') as outfile:
    for fname in filenames:
        if len(set(Path(fname).stem.upper().split('_')) & set(seg_to_keep)) > 0:
            print('yes')
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            outfile.write('\n')

exit(0)
