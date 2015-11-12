# get_BLAT_reads
Scripts to cleanup NGS reads, BLAT to reference and pull the matching reads

Written by Matt Gitzendanner

University of Florida

magitz@ufl.edu

# get_reads_w_blat_hit_pe.py
This script goes through BLAT psl files (without headers, though probably will work if they're there too). It takes the query names that have blat hits, and puts those into a list. It then goes to the original read files and pulls the reads with those names. As it's written now, it assumes there are R1 and R2 files for forward and reverse reads, and a singletons file for unpaired reads.

# run_cut_trim_blat.qsub
This script runs the series of steps to go from raw reads to reads that match a reference.
The steps are:
* Use cutadapt to remove adapters
* Use sickle to filter reads
* Use blat to identify reads with hit to reference(s)
* Use script to pull the reads with hits into new files.
