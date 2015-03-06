#!/bin/perl

use strict ;

my $i ;
my @readFileList ;
my $numOfThread = 1 ;
my $kmerSize = 23 ;
my $bloomFilterSize = 100000000 ;

my $usage = "Usage: perl ./run.pl [OPTIONS]\n".
		"OPTIONS:\n".
		"Required parameters:\n".
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n".
		"\t-k kmer_length\n".
		"Other parameters:\n".
		"\t-od output_file_directory (default: ./)\n".
		"\t-t number of threads to use (default: 1)\n".
		"\t-trim allow trimming (default: false)\n".
		"\t-all: output all the reads including those unfixable (default: false)\n".
		"\t-maxcor: the maximum number of correction every 100bp (default: 8)\n".
		"\t-ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)"; 

if ( scalar( @ARGV ) == 0 )
{
	die "$usage\n" ;
}

for ( $i = 0 ; $i < scalar(@ARGV) ; ++$i )
{
	if ( $ARGV[$i] eq "-r" )
	{
		push @readFileList, $ARGV[$i+1] ;
	}
	elsif ( $ARGV[ $i ] eq "-p" )
	{
		push @readFileList, $ARGV[$i + 1] ;
		push @readFileList, $ARGV[$i + 2] ;
	}
	elsif ( $ARGV[$i] eq "-t" )
	{
		$numOfThread = $ARGV[$i+1] ; 
	}
	elsif ( $ARGV[$i] eq "-k" )
	{
		$kmerSize = $ARGV[$i + 1] ;
	}
	elsif ( $ARGV[$i] eq "-ek" )
	{
		$bloomFilterSize = $ARGV[$i + 1] ;
		# Knock out this argument.
		my $j ; 
		for ( $j = $i ; $j < scalar( @ARGV ) - 2 ; ++$j )
		{
			$ARGV[$j] = $ARGV[$j + 2] ;
		}
		pop @ARGV ;
		pop @ARGV ;
	}
}
#`echo $numOfThread > tmp.out `

print( "Count the kmers\n" ) ;
system( "jellyfish bc -m $kmerSize -s $bloomFilterSize -C -t $numOfThread -o tmp.bc @readFileList" ) ;
system( "jellyfish count -m $kmerSize -s 100000 -C -t $numOfThread --bc tmp.bc -o tmp.mer_counts @readFileList" ) ;
print( "Dump the kmers\n" ) ;
system( "jellyfish dump -L 2 tmp.mer_counts > tmp.jf_dump" ) ;
print( "Error correction\n" ) ;
system( "./rcorrector @ARGV -c tmp.jf_dump" ) ;

#system( "rm tmp.bc tmp.mer_counts tmp.jf_dump" );
