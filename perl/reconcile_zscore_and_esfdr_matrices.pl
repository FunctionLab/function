#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

#################################################################
# This script takes output files from enrichment analyses like  #
# gsea_two-sided_PAGE and then filter for significantly         #
# enriched genesets.                                            #
#                                                               #
# The input files to this script are two matrices, each with    #
# genesets along the rows and conditions along the columns. In  #
# each cell (i,j) the first matrix contains the enrichment      #
# z-scores, and the second contains an enrichemnt score derived #
# from multiple hypothesis corrected q-values, termed ESFDR,    #
# which are typically equal to the signed (-1)*log10(FDR).      #
#                                                               #
# So, based on a user-provided ESFDR cutoff, genesets that pass #
# this cutoff in at least one condition are selected and their  #
# z-scores are reported.                                        #
#################################################################

my($file1, $file2, $esfdrcutoff) = @ARGV;
open FH, "$file1" or die "Can't open $file1!"; chomp(my @f=<FH>); close FH;
open GH, "$file2" or die "Can't open $file2!"; chomp(my @g=<GH>); close GH;
my $outfile = $file1; $outfile =~ s/\.zscore\.mat/\.fdr-based\.zscore\.mat/g;
open HH, ">$outfile";

# File1: Two-sided PAGE Z-score matrix
# File2: Two-sided PAGE ESFDR matrix
# ESFDR cutoff

my $idx_begin = 3; my %esfdr = (); my @p;
for(my $i=1; $i<=$#g; $i++) {
	@p=(); @p=split("\t",$g[$i]);

	for(my $j=$idx_begin; $j<=$#p; $j++) {
		push(@{$esfdr{$p[0]}}, $p[$j]);
	}
}

print HH "$f[0]\n";
for(my $i=1; $i<=$#f; $i++) {
	@p=(); @p=split("\t",$f[$i]);

	my $row_incl_par=0;
	for(my $j=$idx_begin; $j<=$#p; $j++) {
		if(abs(${$esfdr{$p[0]}}[$j-$idx_begin]) >= $esfdrcutoff) {
			$row_incl_par=1; #if($p[$j]<0) { $row_incl_par=-1; }
			last;
		}
		#print HH "\t$row_incl_par";
	}
	#print HH "\n";

	if($row_incl_par eq 1) { print HH "$f[$i]\n"; }
}

close HH;

