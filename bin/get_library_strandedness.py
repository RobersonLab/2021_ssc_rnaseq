#!/usr/bin/env python

###########
# imports #
###########
import argparse
import sys

################################
# setup command-line variables #
################################
parser = argparse.ArgumentParser()

parser.add_argument( 'file', help="path to file" )
parser.add_argument( 'id', help="ID to match" )
parser.add_argument( 'info', help="what to return? number or name", choices=[ 'name', 'number' ] )
parser.add_argument( 'program', help='which program is this for?', choices=[ 'featurecounts', 'salmon' ] )
parser.add_argument( "idtype", help="is this rgid or rgsm?", choices=[ "rgid", "rgsm" ] )

args = parser.parse_args()

#######
# run #
#######
column_dict = {}

response_dict = { 'featurecounts':{'stranded':{'name':'stranded', 'number':1}, 'unstranded':{'name':'unstranded', 'number':0}, 'reverse':{'name':'reverse', 'number':2}}, 'salmon':{'stranded':{'name':"ISF"}, 'unstranded':{'name':"IU"}, 'reverse':{'name':"ISR"}} }

stranded_type_string = None

with open( args.file, 'r' ) as FILE:
	line_num = 0
	
	for line in FILE:
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
			
			file_strandedness = vals[ column_dict[ "Strandedness" ] ]
			
			if args.idtype == "rgid":
				if vals[ column_dict[ 'RGID' ] ] == args.id:
					stranded_type_string = file_strandedness.lower()
					break
			elif args.idtype == "rgsm":
				if vals[ column_dict[ 'RGSM' ] ] == args.id:
					stranded_type_string = file_strandedness.lower()
					break
	
	if stranded_type_string is None:
		sys.exit( 1 )
	
	try:
		print( response_dict[ args.program ][ stranded_type_string ][ args.info ] )
	except:
		sys.exit( 1 )
