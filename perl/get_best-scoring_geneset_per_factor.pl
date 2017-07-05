#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my($in_tab, $in_scol, $out_file) = @ARGV;

# SigPairs table from geneset_enrichment_analysis_of_factored_genelist.pl
# Score to use; Column index; begins at 0
# Outfile: One geneset per factor that has the best score

open FH, "$in_tab" or die "Can't open $in_tab!"; chomp(my @f=<FH>); close FH;
open HH, ">$out_file";

my %gs_max_score = (); my %gs_max_match = (); my @p;
foreach (@f) {
    if($_ =~ /^#/) { print HH "$_\n"; next; }
    @p = split '\t', $_;
    if(exists $gs_max_score{$p[0]}) {
        if($p[$in_scol] > $gs_max_score{$p[0]}) {
            $gs_max_score{$p[0]} = $p[$in_scol];
            $gs_max_match{$p[0]} = $_;
        }
    }
    else {
        $gs_max_score{$p[0]} = $p[$in_scol];
        $gs_max_match{$p[0]} = $_;
    }
}

foreach (sort {$gs_max_score{$b} <=> $gs_max_score{$a}} keys %gs_max_score) {
    print HH "$gs_max_match{$_}\n";
}

close HH;
