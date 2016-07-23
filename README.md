Rcorrector
=========

Described in:

Song, L., Florea, L., [Rcorrector: Efficient and accurate error correction for Illumina RNA-seq reads](http://www.gigasciencejournal.com/content/4/1/48). GigaScience. 2015, 4:48.

Copyright (C) 2012-2013, and GNU GPL, by Li Song and Liliana Florea

Rcorrrector includes the program [Jellyfish2](http://www.genome.umd.edu/jellyfish.html)

### What is Rcorrector?

Rcorrector(RNA-seq error CORRECTOR) is a kmer-based error correction method for RNA-seq data. 

Rcorrector can also be applied to other type of sequencing data where the read coverage is non-uniform, such as single-cell sequencing.

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/rcorrector), e.g. with `git clone https://github.com/mourisl/rcorrector.git`
2. Run `make` in the repo directory
	During the `make` procedure, the script will check whether you have jellyfish2 in $PATH. If not, it will download and compile jellyfish2 from its repository. 

### Usage
	Usage: perl run_rcorrector.pl [OPTIONS]
	OPTIONS:
		Required
		-s seq_files: comma separated files for single-end data sets
		-1 seq_files_left: comma separated files for the first mate in the paried-end data sets
		-2 seq_files_right: comma separated files for the second mate in the paired-end data sets
		-i seq_files_interleaved: comma sperated files for interleaved paired-end data sets
		Optional
		-k INT: kmer_length (<=32, default: 23)
		-od STRING: output_file_directory (default: ./)
		-t INT: number of threads to use (default: 1)
		-trim : allow trimming (default: false)
		-maxcorK INT: the maximum number of correction within k-bp window (default: 4)
		-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold, lower for more divergent genome (default: 0.95)
		-ek INT: expected number of kmers; does not affect the correctness of program but affects the memory usage (default: 100000000)
		-stdout: output the corrected reads to stdout (default: not used)
		-verbose: output some correction information to stdout (default: not used)
		-stage INT: start from which stage (default: 0)
			0-start from begining(storing kmers in bloom filter) ;
			1-start from count kmers showed up in bloom filter;
			2-start from dumping kmer counts into a jf_dump file;
			3-start from error correction.


### Output
For each input file, Rcorrector will generate the corresponding output file with "*.cor.fq/fa" in the directory specified by "-od". 

In the header line for each read, Rcorrector will append some information.

	"cor": some bases of the sequence are corrected
	"unfixable_error": the errors could not be corrected
	"l:INT m:INT h:INT": the lowest, median and highest kmer count of the kmers from the read


### Example
We put a small sample data set, you can run them by:

perl run_rcorrector.pl -1 Sample/sample_read1.fq -2 Sample/sample_read2.fq  

### Contact
lsong10@jhu.edu

### Terms of use

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General
Public License along with this program; if not, you can obtain one from
http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
### Support

Create a [GitHub issue](https://github.com/mourisl/rcorrector/issues).


