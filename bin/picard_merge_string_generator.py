#!/usr/bin/env python

import argparse
import os
import sys

parser = argparse.ArgumentParser()

parser.add_argument( 'file', help="path to file" )
parser.add_argument( 'RGSM', help="RGSM for sample to merge" )
parser.add_argument( '--starpath', help="Both path to input BAMs and path for output", default="bam/" )
parser.add_argument( 'mode', choices=['mergecmd', 'rmcmd' ] )

args = parser.parse_args()

buffer = ""
column_dict = {}

sample_list = []

file_line = 0

with open( args.file, 'r' ) as FILE:
	for line in FILE:
		file_line += 1
		line = line.rstrip()
	
		if file_line == 1:
			line = line.lstrip( '#' )
			
			vals = line.split( ',' )
			
			col_index = -1
			for curr_val in vals:
				col_index += 1
				column_dict[ curr_val ] = col_index
		else:
			vals = line.split( ',' )
			
			if vals[ column_dict[ 'RGSM' ] ] == args.RGSM:
				sample_list.append( vals[ column_dict[ 'RGID' ] ] )

sample_list = [ os.path.normpath( "%s/%s/Aligned.out.bam" % ( args.starpath, s ) ) for s in sample_list ]
				
try:
	inputs = 'INPUT=%s' % ( ' INPUT='.join( sample_list ) )
	removes = ' '.join( sample_list )
	
	if args.mode == 'mergecmd':
		#print( "%s" % ( inputs ), end="" )
		print( inputs, end="" )
	elif args.mode == 'rmcmd':
		#print( "%s" % ( removes ), end="" )
		print( removes, end="" )
	else:
		sys.exit( 1 )
except:
	sys.exit( 1 )

