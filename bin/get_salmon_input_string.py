#!/usr/bin/env python

###########
# imports #
###########
import argparse
import os
import sys

################################
# setup command-line variables #
################################
parser = argparse.ArgumentParser()

parser.add_argument( 'rgsm', help='path to file of rgids, one per line' )
parser.add_argument( 'rginfo_file', help="path to rginfo file" )
parser.add_argument( '--fastqpath', help="What is the relative or absolute path to all FASTQ files?", default="." )
parser.add_argument( '--r1_ext', help="What to append to the end of all R1s", default="_1_clean.fastq.gz", required=False )
parser.add_argument( '--r2_ext', help="What to append to the end of all R2s", default="_2_clean.fastq.gz", required=False )

args = parser.parse_args()

args.fastqpath = os.path.normpath( args.fastqpath )
		
################################
# process read group info file #
################################
rgid_list = []

with open( args.rginfo_file, 'r' ) as RGINFO:
	column_dict = {}
	line_num = 0
	
	for line in RGINFO:
		line = line.rstrip()
		
		if len( line ) == 0:
			continue
			
		line_num += 1
		
		if line_num == 1:
			line = line.lstrip( '#' )
			
			vals = line.split( ',' )
			
			col_index = -1
			for curr_val in vals:
				col_index += 1
				column_dict[ curr_val ] = col_index
		else:
			vals = line.split( ',' )
			curr_rgsm = vals[ column_dict[ 'RGSM' ] ]
			curr_rgid = vals[ column_dict[ 'RGID' ] ]
			
			if curr_rgsm == args.rgsm:
				rgid_list.append( curr_rgid )

rgid_list.sort()

r1_list = []
r2_list = []

for rgid in rgid_list:
	r1_list.append( os.path.join( args.fastqpath, "%s%s" % ( rgid, args.r1_ext ) ) )
	r2_list.append( os.path.join( args.fastqpath, "%s%s" % ( rgid, args.r2_ext ) ) )
	
out_string = "-1 %s -2 %s" % ( " ".join( r1_list ), " ".join( r2_list ) )

print( out_string )
