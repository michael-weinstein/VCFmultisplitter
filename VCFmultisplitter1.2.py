#!/usr/bin/env python

#This program will create individual VCFs from a multi-patient VCF.
#Version information for 1.1: This version introduces an option to split for de novo searches.  This option
#will keep all information, including null reads or homozygous reference reads.
#These reads will otherwise not be included in the final outputs.
#Copyright 2014, Michael Weinstein, Daniel Cohn laboratory

def fileinfo(file):  #subroutine designed to get information on the VCF to be processed
    import re  #import the library needed to use regex
    vcf = open (file, 'r') #opens the VCF to be processed
    headerlines = '' #initializes headerlines as a null string
    currentline = vcf.readline()  #reads a line of the VCF into the currentline string
    currentline = currentline.strip('"')  #removes any leading or trailing quotation marks that may have been added to the file by a third-party program used to view the VCF (MS Excel, I'm looking in your direction here)
    while currentline:  #so long as there is some value in the currentline string (meaning it was not a blank file and we have not hit the end of the file)
        if re.match ('^\#\#.*$', currentline): #uses a regex to see if the currentline starts with two hash tags
            headerlines += currentline  #if so, it gets added to the growing string of headerlines
        elif re.match('^\#CHROM.*$', currentline):  #if not, it should start with #CHROM, indicating that it is the column header row
            columnheaders = currentline  #if so, it reads the line into columnheaders
            columns = currentline.split ('\t')  #splits the columnheaders string on tabs to create an array called columns
            genomes = len(columns) - 9  #this assumes the standard 9 columns that should be found in every VCF and subtracts it from the total number of columns to count the number of individuals on the VCF
            vcf.close()  #closes the VCF file
            return (genomes, columnheaders, headerlines)  #returns a tupple containing the number of genomes, the column headers, and the headerlines while exiting this subroutine
        else:  #if the line does not match either of the above line types, something is very wrong with the VCF
            usage('VCF malformation in the header lines.')  #prints an error message and instructions
            quit() #and drops out of the program
        currentline = vcf.readline()  #reads the next line into currentline before reiterating through the loop
        currentline = currentline.strip('"')  #removes any leading or trailing quotes that may have been added to the line by a third-party program used to view the VCF
    usage('VCF malformation in or near headerlines, possible missing data.')  #if we exit the loop before we exit the subroutine, something is very wrong (we hit the end of file before finding all the header rows and column headers)
    quit()  #so we print an error message and quit
    
def checkargs():  #subroutine for validating commandline arguments
    import argparse #loads the required library for reading the commandline
    import os  #imports the library we will need to check if the file exists
    parser = argparse.ArgumentParser()
    parser.add_argument ("-f", "--file", help = "Specify the desired file to cut for submission.")  #tells the parser to look for -f and stuff after it and call that the filename
    parser.add_argument ("-d", "--denovo", help = "De novo mode, keep variants that are present in other VCFs, but not the one being listed", action = 'store_true')
    args = parser.parse_args()  #puts the arguments into the args object
    if not args.file:  #if the args.file value is null, give an error message and quit the program
        usage("No file specified.") 
        quit()
    elif not os.path.isfile(args.file):  #if the file specified in the arguments doesn't exist, quit the program and give an error message
        usage("Could not locate " + args.file + "on this system.") 
        quit()
    else:
        return (args.file,args.denovo)  #returns the validated filename to the main program
    
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
denovos = False    
args = checkargs()
filename = args[0]  #puts the command line arguments into the filename string after validating it (see checkargs subroutine)
denovos = args[1]
fileinfotupple = fileinfo(filename)  #number of cases on the vcf, the column titles, and the header (info) lines
ncases = fileinfotupple[0]  #This and the next two lines assign values from the returned tupple of file info to properly-named variables
infolines = fileinfotupple[2]
columnnames = fileinfotupple[1].split('\t') #splits the columnnames on tabs into an array as it recovers it from the file info tupple
hardcolumns = []  #intializes an empty array to add column header names to in order
for i in range (0,9):  #iterates through the standard VCF column names
    hardcolumns.insert(i,columnnames[i])  #and adds them to the list of hard columns
casearray = [] #initializes an empty list for case names on the VCF
for i in range (9, len(columnnames)): #iterates through the column header list starting at position 9 where the case names begin through the last case name
    casearray.insert(i-9,columnnames[i].rstrip()) #and reads them into an array while removing any trailing end of line characters that show up
outputfilearray = [] #initializes an empty array for the list of output files
for i in range (0, ncases): #iterates through the new list of cases
    outputfilearray.insert(i,filename + '.split.' + casearray[i] + '.vcf')  #creates an array of filenames based on the cases in the VCF, one file per case
if not filenamesfree(outputfilearray):  #runs a subroutine to check if the filenames that might be needed are free (avoids accidental overwrites of existing data)
    usage('Potential output file name may already be taken.') #if the filenames are already taken, returns an error message and instructions
    quit()
for i in range (0, ncases): #iterates through a loop based on the number of cases
    output = open(outputfilearray[i],'w') #creates and opens an output file for each case
    output.write (infolines) #and writes the infolines (the first set of headerlines containing the VCF information) to each file
    for item in hardcolumns:    #iterates through the standard VCF column array
        output.write (item + '\t') #and writes each one to the appropriate position in the file to generate column headers
    output.write (casearray[i] + '\n') #writes the name of the case to the final spot on the column header line
    output.close()  #closes file and returns to the beginning of the loop
vcf = open(filename,'r')  #opens the vcf file we will be processing
progress = 0  #initializes our progress counter to 0
answer = False  #initializes our variable to test for a legitimate answer to False
while not answer and not denovos:  #enters the loop and stays in it until answer is equal to True
    print ('Will this VCF be used for finding potential de novo variants? (Y/N)')
    answer = input('>>') #sets answer equal to some value input by the user
    if str(answer) == 'y' or str(answer) == 'Y':  #checks if the answer is a yes
        denovos = True #if so, sets denovos to True (tells the program we are looking for denovos)
    elif str(answer) == 'n' or str(answer) == 'N': #if the answer is No
        denovos = False  #sets the denovos mode to False (neither of these options change answer to false, so the loop can exit)
    else: #if the answer is not a value indicating a yes or no
        print ('Invalid response.')
        answer = False #set ansewr to false so the loop will continue until a satisfactory answer is given
line = vcf.readline()  #this is a throw-away line, as we already have the header string stored, but at least initializes line to a non-null value
while line:  #checks to make sure that there is a value in line (and we are not at end of file)
    progress +=1 #increment the progress counter
    print ('Processing line ' + str(progress) + '.\r', end = '') #display the progress counter to the user and update it
    line = vcf.readline() #read a new line from the file
    if not line: #if the line is null (indicating end of file)
        continue #return to the beginning of the loop, the null value will make it exit
    line = line.strip('"') #remove leading and trailing quotes from the line
    if re.match('^#.*$', line): #check if the line starts with a hastag
        continue #if so, returns to the beginning of the loop to read the next line and test it
    linearray = line.split('\t') #if the line is not one of the headerlines (as checked for above), this will split the line into an array on the tabs
    for i in range (0, len(linearray)): #this will iterate through each value in the line array
        linearray[i] = linearray[i].strip() #this will remove any leading or trailing spaces or end of lines
        linearray[i] = linearray[i].strip('"')  #this will remove any leading or trailing quotes
    for i in range (0,ncases):  #this will start a loop to iterate through each case
        regex = re.search('^(./.).*$', linearray[i+9]) #searches the individual we are looking at based on our position in the loop and finds (any character/any character) which would indicate their genotype
        genotype = regex.group(1)  #uses the value returned by the regex to determine their genotype
        if not denovos and (genotype == '0/0' or genotype == './.'):  #first checks to see if we are looking for de novo variants (as indicated by the user) with this VCF, then checks to see if the read was either a null read or a homozygous reference read
            continue  #if we are not looking for de novo variants and the read was either reference or a null read, we continue the loop without writing (as there is no variant to report here)
        output = open(outputfilearray[i],'a')  #otherwise, we open the appropriate output file
        for j in range (0,9):  #and iterate through the standard VCF columns
            output.write (linearray[j] + '\t') #and write each column to the patient's file
        output.write (linearray[i+9] + '\n') #followed by the individual's read data
        output.close()  #then close the output file and move on to the next
vcf.close() #this closes the input VCF after everything is done
print ('\nDone!\n') #and reports to the user that we are done before quitting