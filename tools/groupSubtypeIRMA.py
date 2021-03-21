import os
import sys
import subprocess
import glob

# python groupSubtype.py assemblies/{sample}/contigs.fasta workspace
# gotta recode this piece of shit. -- ammar

filename = sys.argv[1]
with open(filename, 'r') as fi:
    fluType = ""
    subtypeHA = ""
    subtypeNA = ""
    for line in fi:
        if ">B" in line:
            fluType = "B"
        elif ">A" in line:
            fluType = "A"
            lineclean = line.upper()
            if "_HA_" in lineclean or "_NA_" in lineclean:
                if "_HA_" in lineclean:
                    subtypeHA = line.split("_")[2][:-1]
                elif "_NA_" in lineclean:
                    subtypeNA = line.split("_")[2][:-1]
    if fluType == "B":
        bname = filename.split("/")[-2]
        # Create folder
        command = "mkdir -p " + sys.argv[2] + "assemblies/rename/type/FLUB/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        # Move file
        command = "cp " + sys.argv[2] + "assemblies/rename/"+bname + \
            ".fa " + " " + sys.argv[2] + "assemblies/rename/type/FLUB/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
    elif fluType == "A":
        bname = filename.split("/")[-2]
        # Create folder
        command = "mkdir -p " + sys.argv[2] + "assemblies/rename/type/FLUA/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        # Create folder
        command = "mkdir " + sys.argv[2] + "assemblies/rename/type/FLUA/"+subtypeHA+subtypeNA+"/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        # Move file
        command = "cp " + sys.argv[2] + "assemblies/rename/"+bname+".fa " + \
            sys.argv[2] + "assemblies/rename/type/FLUA/"+subtypeHA+subtypeNA+"/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
