#!/usr/bin/env python

# syntax written with Python 3 in mind

##########
# Import #
##########
import argparse
from builtins import int
from builtins import str
from functools import total_ordering
import gzip
import sys

####################
# Version and name #
####################
__script_path__ = sys.argv[0]
__script_name__ = __script_path__.split( '/' )[-1].split( '\\' )[-1]
__version__ = '1.0.0'

################
# SJ tab class #
################
@total_ordering
class RNAstarSpliceJunction:
	"""A class to deal with RNA-STAR splice junction tab data"""
	__slots__ = [ 'chrom', 'intronStart', 'intronEnd', 'strand', 'motif', 'annotationStatus', 'uniqueReads', 'multimapReads', 'maxSpliceOverhang' ]
	
	def __init__( self, value_array ):
		# 0 chrom
		# 1 intron starts
		# 2 intron ends
		# 3 strand
		# 4 motif
		# 5 annotation status
		# 6 unique reads crossing
		# 7 multipmap reads crossing
		# 8 maximum splice alignment overhang
		
		self.chrom = value_array[ 0 ]
		self.intronStart = int( value_array[ 1 ] )
		self.intronEnd = int( value_array[ 2 ] )
		self.strand = value_array[ 3 ]
		self.motif = int( value_array[ 4 ] )
		self.annotationStatus = int( value_array[ 5 ] )
		self.uniqueReads = int( value_array[ 6 ] )
		self.multimapReads = int( value_array[ 7 ] )
		self.maxSpliceOverhang = int( value_array[ 8 ] )
		
		################
		# setup strand #
		################
		if self.strand == '1':
			self.strand = '+'
		elif self.strand == '2':
			self.strand = '-'
		elif self.strand == '0':
			self.strand = '.'
		
		#####################
		# validity checking #
		#####################
		if self.intronStart < 0:
			print( "Tried to set an intronStart < 0 [%s]" % self.intronStart )
			sys.exit( 1 )
			
		if self.intronEnd < 0:
			print( "Tried to set an intronEnd < 0 [%s]" % self.intronEnd )
			sys.exit( 1 )
		
		assert( self.motif in ( 0, 1, 2, 3, 4, 5, 6 ) )
		assert( self.annotationStatus in ( 0, 1 ) )
	
	#########################
	# other class functions #
	#########################
	def is_new_jxn( self ):
		if self.annotationStatus == 0:
			return True
		return False
	
	def set_counts( self, value=None ):
		if value is None:
			value = 1
		
		self.uniqueReads = value
		self.multimapReads = value
	
	def __len__( self ):
		return self.intronEnd - self.intronStart + 1
		
	def __add__( self, other ):
		unique_count = self.uniqueReads + other.uniqueReads
		multi_count = self.multimapReads + other.multimapReads
		
		newMaxSpliceOverhang = max( self.maxSpliceOverhang, other.maxSpliceOverhang )
		
		return RNAstarSpliceJunction( [ self.chrom, self.intronStart, self.intronEnd, self.strand, self.motif, self.annotationStatus, unique_count, multi_count, newMaxSpliceOverhang ] )
	
	def __str__( self ):
		return "%s_%s_%s_%s" % ( self.chrom, self.intronStart, self.intronEnd, self.strand )
		
	def __repr__( self ):
		return str( self )
		
	def print_tab( self ):
		return "%s\t%s\t%s\t%s" % ( self.chrom, self.intronStart, self.intronEnd, self.strand )
		#return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ( self.chrom, self.intronStart, self.intronEnd, self.strand, self.motif, self.annotationStatus, self.uniqueReads, self.multimapReads, self.maxSpliceOverhang )
		
	def __eq__( self, other ):
		return ( self.chrom == other.chrom and self.intronStart == other.intronStart and self.intronEnd == other.intronEnd )
		
	def __lt__( self, other ):
		if self.chrom < other.chrom:
			return True
		elif self.chrom == other.chrom and self.intronStart < other.intronStart:
			return True
		elif self.chrom == other.chrom and self.intronStart == other.intronStart and self.intronEnd < other.intronEnd:
			return True
		else:
			return False
			
	def lt( self, other ):
		self.__lt__( other )
		
	def eq( self, other ):
		self.__eq__( other )
		
#####################
# main run function #
#####################
def run():
	#######################
	# Set argparse values #
	#######################
	parser = argparse.ArgumentParser( prog=__script_name__, epilog="%s v%s" % ( __script_name__, __version__ ) )
	
	parser.add_argument( "fastaindex", help="Path to FASTA index file corresponding to the genome use for RNASTAR alignment" )
	parser.add_argument( "sjtab_filelist", help="Path to file containing paths (relative to script execution or absolute) to the SJ tab files to merge" )
	parser.add_argument( "--mincount", help="only print sites with at least this many reads", default=1, type=int )
	parser.add_argument( "--track_known", help="Whether to add known junctions to the dictionary", default=False, action="store_true" )

	input_args = parser.parse_args()

	#################
	# List of files #
	#################
	file_list = []
	
	with open( input_args.sjtab_filelist, 'r' ) as SJLIST:
		for line in SJLIST:
			line = line.rstrip( "\r\n" )
			
			if len( line ) > 0:
				file_list.append( line )
			
	file_list = list( set( file_list ) )
	
	###############################
	# Process FAI object in order #
	###############################
	contig_list = []
	
	with open( input_args.fastaindex, 'r' ) as FAI:
		for line in FAI:
			line = line.rstrip( "\r\n" )
			
			values = line.split( "\t" )
			
			contig_list.append( values[0] )
			
	contig_list.sort()
	
	#################
	# Wrangle files #
	#################
	splicejxn_dict = {}
	
	for curr_sjtabfile in file_list:
		with gzip.open( curr_sjtabfile, 'rt' ) as SJTAB:
			for line in SJTAB:
				line = line.strip( "\r\n" )
				
				curr_sj_obj = RNAstarSpliceJunction( line.split( '\t' ) )
				
				if not input_args.track_known and not curr_sj_obj.is_new_jxn():
					continue
				
				sj_name = str( curr_sj_obj )
				sj_chrom = curr_sj_obj.chrom
				
				if sj_chrom in splicejxn_dict:
					if sj_name in splicejxn_dict[ sj_chrom ]:
						splicejxn_dict[ sj_chrom ][ sj_name ] = splicejxn_dict[sj_chrom][ sj_name ] + curr_sj_obj
					else:
						splicejxn_dict[ sj_chrom ][ sj_name ] = curr_sj_obj
				else:
					splicejxn_dict[ sj_chrom ] = { sj_name:curr_sj_obj }
					
	##################################################################
	# By chromosome, in order, convert to values and print to screen #
	##################################################################
	for curr_contig in contig_list:
		if curr_contig in splicejxn_dict:
			contig_sj_dict = splicejxn_dict[ curr_contig ]
			
			for curr_junction in sorted( contig_sj_dict.values() ):
				if not input_args.track_known:
					if not curr_junction.is_new_jxn():
						continue
					
				if curr_junction.uniqueReads >= input_args.mincount:
					print( curr_junction.print_tab() )

if __name__ == "__main__":
	run()
