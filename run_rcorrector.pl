#!/bin/perl

use strict ;
use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename;

my $i ;
my @readFileList ;
my $numOfThread = 1 ;
my $kmerSize = 23 ;
my $bloomFilterSize = 100000000 ;
my $WD = dirname( abs_path( $0 ) ) ;


my $usage = "Usage: perl ./run_rcorrector.pl [OPTIONS]\n".
		"OPTIONS:\n".
		"Required parameters:\n".
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n".
		"\t-p seq_file_left seq_file_right: the paths to the paired-end data set. Can use multiple -p to specifiy multiple sequence files\n".
		"Other parameters:\n".
		"\t-k kmer_length (default: 23)\n".
		"\t-od output_file_directory (default: ./)\n".
		"\t-t number of threads to use (default: 1)\n".
		"\t-trim allow trimming (default: false)\n".
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

system( "$WD/jellyfish-2.1.3/bin/jellyfish bc -m $kmerSize -s $bloomFilterSize -C -t $numOfThread -o tmp.bc @readFileList" ) ;
system( "$WD/jellyfish-2.1.3/bin/jellyfish count -m $kmerSize -s 100000 -C -t $numOfThread --bc tmp.bc -o tmp.mer_counts @readFileList" ) ;
print( "Dump the kmers\n" ) ;
system( "$WD/jellyfish-2.1.3/bin/jellyfish dump -L 2 tmp.mer_counts > tmp.jf_dump" ) ;
print( "Error correction\n" ) ;
system( "$WD/rcorrector @ARGV -c tmp.jf_dump -k $kmerSize" ) ;

#system( "rm tmp.bc tmp.mer_counts tmp.jf_dump" );
