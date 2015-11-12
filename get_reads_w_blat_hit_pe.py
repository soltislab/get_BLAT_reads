#!/usr/bin/env python
import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
import os, errno

# =====================================================
#  
# This script goes through BLAT psl files (without headers, though probably will work if they're there too)
# It takes the query names that have blat hits, and puts those into a list.
# It then goes to the original read files and pulls the reads with those names.
# As it's written now, it assumes there are R1 and R2 files for forward and reverse reads, and a singletons 
# file for unpaired reads.
#
#  Matt Gitzendanner
#  University of Florida
#  magitz@ufl.edu
#  2/19/14
#
# =====================================================


#####################
# Parse commandline options.
#  -b blat filename name-- without the R1.psl, R2.psl, singletons.psl
#  -r Directory with reads files-- R1, R2, singletons
#  -o Output path
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-b", help="blat name-- without the R1.psl, R2.psl, singletons.psl")
parser.add_argument("-r", help="Directory with reads files-- R1, R2, singletons")
parser.add_argument("-o", help="Output path")
args = parser.parse_args()

blat_name= args.b
reads_dir=args.r
out_dir=args.o


#Function to go through a psl file and add the hits to an array.
# Input: psl file path
def add_to_list(psl_file):
	try:
		BLAT=open(psl_file, 'r')
	except IOError:
		print "Can't open file:", psl_file

	for Line in BLAT :
		Line = Line.strip('\n')
		Line_bits=re.split('\t',Line)
	
		read_name=Line_bits[9]

		if not (read_name in read_list):	#add the scaffold to the list for the taxon.
			read_list.append(read_name)

#Function to go through a reads file and pull reads that are in a list.
#Input: Reads file path and output file path			
def get_reads(reads_file,out_file):
	try:
		READS=open(reads_file, 'r')
	except IOError:
		print "Can't open file:", reads_file

	try:
		OUT=open(out_file, 'w')
	except IOError:
		print "Can't open file:", out_file
	
	for record in SeqIO.parse(READS, "fastq") :
		if record.id in read_list:
			SeqIO.write(record, OUT, "fastq")
			
#Make the output directory if it's not already made.
try:
	os.makedirs(out_dir)
except OSError, err:
    # Reraise the error unless it's about an already existing directory 
    if err.errno != errno.EEXIST or not os.path.isdir(out_dir): 
        raise


read_list=[]	#List for the reads that need to be pulled


#Run add_to_list on the 3 psl files.
blat_file=blat_name + "R1.psl"
add_to_list(blat_file)

blat_file=blat_name + "R2.psl"
add_to_list(blat_file)

blat_file=blat_name + "singletons.psl"
add_to_list(blat_file)



#Run get_reads on the 3 original reads files.
file=os.path.basename(blat_name)
	
reads_file= reads_dir + "/" + file + "R1.fq"
out_file= out_dir + "/" + file + "R1.fq"
get_reads(reads_file, out_file)

reads_file= reads_dir + "/" + file + "R2.fq"
out_file= out_dir + "/" + file + "R2.fq"
get_reads(reads_file, out_file)

reads_file= reads_dir + "/" + file + "singletons.fq"
out_file= out_dir + "/" + file + "singletons.fq"
get_reads(reads_file, out_file)
