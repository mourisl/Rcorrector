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
		"\t-s seq_files: comma separated files for single-end data sets\n".
		"\t-1 seq_files_left: comma separated files for the first mate in the paried-end data sets\n".
		"\t-2 seq_files_right: comma separated files for the second mate in the paired-end data sets\n".
		"Other parameters:\n".
		"\t-k kmer_length (default: 23)\n".
		"\t-od output_file_directory (default: ./)\n".
		"\t-t number of threads to use (default: 1)\n".
		#"\t-trim allow trimming (default: false)\n".
		"\t-maxcor: the maximum number of correction every 100bp (default: 8)\n".
		"\t-ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)"; 

if ( scalar( @ARGV ) == 0 )
{
	die "$usage\n" ;
}

my $fileArguments ;
my @singleFileList ;
my @firstMateFileList ;
my @secondMateFileList ;

my @rcorrectorArguments ;

for ( $i = 0 ; $i < scalar(@ARGV) ; ++$i )
{
	if ( $ARGV[$i] eq "-r" )
	{
		push @readFileList, $ARGV[$i+1] ;

		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[ $i ] eq "-p" )
	{
		push @readFileList, $ARGV[$i + 1] ;
		push @readFileList, $ARGV[$i + 2] ;
		
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;
		push @rcorrectorArguments, $ARGV[$i + 2] ;

		$i += 2 ;
	}
	elsif ( $ARGV[ $i ] eq "-s" )
	{
		my @cols = split /,/, $ARGV[$i + 1] ;
		my $j ;
		for ( $j = 0 ; $j < scalar( @cols ) ; ++$j )
		{
			push @readFileList, $cols[ $j ] ;
			push @singleFileList, $cols[ $j ] ;
		}
		++$i ;
	}
	elsif ( $ARGV[ $i ] eq "-1" )
	{
		my @cols = split /,/, $ARGV[$i + 1] ;
		my $j ;
		for ( $j = 0 ; $j < scalar( @cols ) ; ++$j )
		{
			push @readFileList, $cols[ $j ] ;
			push @firstMateFileList, $cols[ $j ] ;
		}
		++$i ;
	}
	elsif ( $ARGV[ $i ] eq "-2" )
	{
		my @cols = split /,/, $ARGV[$i + 1] ;
		my $j ;
		for ( $j = 0 ; $j < scalar( @cols ) ; ++$j )
		{
			push @readFileList, $cols[ $j ] ;
			push @secondMateFileList, $cols[ $j ] ;
		}
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-t" )
	{
		$numOfThread = $ARGV[$i+1] ; 
		
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-k" )
	{
		$kmerSize = $ARGV[$i + 1] ;
		
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-ek" )
	{
		$bloomFilterSize = $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-maxcor" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-maxcorK" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	else
	{
		die "Unknown argument ".$ARGV[$i]."\n" ;
	}
}
#`echo $numOfThread > tmp.out `

# Build the input file arguments for 
for ( my $i = 0 ; $i < @singleFileList ; ++$i )
{
	$fileArguments = $fileArguments." -r ".$singleFileList[$i] ;
}


die "The number of files from -1,-2 should be the same" if ( scalar( @firstMateFileList ) != scalar( @secondMateFileList ) ) ;
for ( my $i = 0 ; $i < @firstMateFileList ; ++$i )
{
	$fileArguments = $fileArguments." -p ".$firstMateFileList[$i]." ".$secondMateFileList[$i] ;
}

print( "Count the kmers\n" ) ;

print( "$WD/jellyfish-2.1.3/bin/jellyfish bc -m $kmerSize -s $bloomFilterSize -C -t $numOfThread -o tmp.bc @readFileList\n" ) ;
system( "$WD/jellyfish-2.1.3/bin/jellyfish bc -m $kmerSize -s $bloomFilterSize -C -t $numOfThread -o tmp.bc @readFileList" ) ;
print( "$WD/jellyfish-2.1.3/bin/jellyfish count -m $kmerSize -s 100000 -C -t $numOfThread --bc tmp.bc -o tmp.mer_counts @readFileList\n" ) ;
system( "$WD/jellyfish-2.1.3/bin/jellyfish count -m $kmerSize -s 100000 -C -t $numOfThread --bc tmp.bc -o tmp.mer_counts @readFileList" ) ;

print( "Dump the kmers\n" ) ;
print( "$WD/jellyfish-2.1.3/bin/jellyfish dump -L 2 tmp.mer_counts > tmp.jf_dump\n" ) ;
system( "$WD/jellyfish-2.1.3/bin/jellyfish dump -L 2 tmp.mer_counts > tmp.jf_dump" ) ;

print( "Error correction\n" ) ;
print( "$WD/rcorrector @rcorrectorArguments $fileArguments -c tmp.jf_dump\n" ) ;
system( "$WD/rcorrector @rcorrectorArguments $fileArguments -c tmp.jf_dump" ) ;

#system( "rm tmp.bc tmp.mer_counts tmp.jf_dump" );
