#!/usr/bin/env python

#This program will create individual VCFs from a multi-patient VCF.
#Copyright 2014, Michael Weinstein, Daniel Cohn laboratory

def fileinfo(file):
    import re
    vcf = open (file, 'r')
    headerlines = ''
    currentline = vcf.readline()
    currentline = currentline.strip('"')
    while currentline:
        if re.match ('^\#\#.*$', currentline):
            headerlines += currentline
        elif re.match('^\#CHROM.*$', currentline):
            columnheaders = currentline
            columns = currentline.split ('\t')
            genomes = len(columns) - 9
            vcf.close()
            return (genomes, columnheaders, headerlines)
        else:
            usage('VCF malformation in the header lines.')
            quit()
        currentline = vcf.readline()
        currentline = currentline.strip('"')
    
def checkargs():  #subroutine for validating commandline arguments
    import argparse #loads the required library for reading the commandline
    import os  #imports the library we will need to check if the file exists
    parser = argparse.ArgumentParser()
    parser.add_argument ("-f", "--file", help = "Specify the desired file to cut for submission.")  #tells the parser to look for -f and stuff after it and call that the filename
    args = parser.parse_args()  #puts the arguments into the args object
    if not args.file:  #if the args.file value is null, give an error message and quit the program
        usage("No file specified.") 
        quit()
    elif not os.path.isfile(args.file):  #if the file specified in the arguments doesn't exist, quit the program and give an error message
        usage("Could not locate " + args.file + "on this system.") 
        quit()
    else:
        return str(args.file)  #returns the validated filename to the main program
    
def filenamesfree(files):  #this subroutine checks if the series of filenames we are likely to need is free
    import os  #imports the library we will need to check filenames in the directory
    for file in files:  #creates a loop to check for 10 possible filenames we can use
        if os.path.isfile(file):  #checks each potential output filename as the loop iterates
            return False  #if it finds that a file by the current name being tested exists, exits the subroutine returning a false value
    return True  #if it finds none of the files exist, it returns a true value

def usage(sin):  #This subroutine prints directions
    print ('Error:' + sin)
    print ('This script will break a VCF that is too large for upload to SeattleSeq into smaller files that can be submitted.')
    print ('Sample commandline:\npython seattleslicer.py -f file.csv')

import re
    
filename = checkargs()  #puts the command line arguments into the filename string after validating it (see checkargs subroutine)
fileinfotupple = fileinfo(filename)  #number of cases on the vcf, the column titles, and the header (info) lines
ncases = fileinfotupple[0]
infolines = fileinfotupple[2]
columnnames = fileinfotupple[1].split('\t')
hardcolumns = []
for i in range (0,9):
    hardcolumns.insert(i,columnnames[i])
casearray = []
for i in range (9, len(columnnames)):
    casearray.insert(i-9,columnnames[i].rstrip())
outputfilearray = []
for i in range (0, ncases):
    outputfilearray.insert(i,filename + '.split.' + casearray[i] + '.vcf')
if not filenamesfree(outputfilearray):  #runs a subroutine to check if the filenames that might be needed are free (avoids accidental overwrites of existing data)
    usage('Potential output file name may already be taken.') #if the filenames are already taken, returns an error message and instructions
    quit()
for i in range (0, ncases):
    output = open(outputfilearray[i],'w')
    output.write (infolines)
    for item in hardcolumns:    
        output.write (item + '\t')
    output.write (casearray[i] + '\n')
    output.close() 
vcf = open(filename,'r')
progress = 0
line = vcf.readline()  #this is a throw-away line, as we already have the header string stored, but at least initializes line to a non-null value
while line:
    progress +=1
    print ('Processing line ' + str(progress) + '.\r', end = '')
    line = vcf.readline()
    if not line:
        continue
    line = line.strip('"')
    if re.match('^#.*$', line):
        continue
    linearray = line.split('\t')
    for i in range (0, len(linearray)):
        linearray[i] = linearray[i].strip()
        linearray[i] = linearray[i].strip('"')
    for i in range (0,ncases):
        regex = re.search('^(./.).*$', linearray[i+9])
        genotype = regex.group(1)
        if genotype == '0/0' or genotype == './.':
            continue
        output = open(outputfilearray[i],'a')
        for j in range (0,9):
            output.write (linearray[j] + '\t')
        output.write (linearray[i+9] + '\n')
        output.close()    
vcf.close()
print ('\nDone!\n')