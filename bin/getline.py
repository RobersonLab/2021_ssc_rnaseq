#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument( 'file', help="path to file" )
parser.add_argument( 'linenumber', help="line to pull", type=int )

args = parser.parse_args()

if args.linenumber < 1:
	print( "Line specified is less than 1. Exiting" )
	sys.exit( 1 )

line_count = 0
buffer = ""

with open( args.file, 'r' ) as FILE:
	for buffer in FILE:
		line_count += 1
		
		if line_count == args.linenumber:
			buffer = buffer.rstrip()
			
			if len( buffer ) == 0:
				print( "Empty line" )
				sys.exit( 1 )
			else:
				print( buffer )
				sys.exit( 0 )

print( "Line beyond end of file" )
sys.exit( 1 )
