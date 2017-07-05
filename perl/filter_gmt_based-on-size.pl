#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my($in_gmt, $min_genes, $max_genes, $out_gmt) = @ARGV;

open GMT, "$in_gmt" or die "Can't open $in_gmt!";
    chomp(my @gmt=<GMT>); close GMT;
open OUT, ">$out_gmt";

my (@p, $ngenes);
my %all_genes = (); my $tot_gs = 0;
foreach (@gmt) {
    @p = split '\t', $_;

    $ngenes = scalar(@p)-2;
    if(($min_genes <= $ngenes)
            and ($ngenes <= $max_genes)) {

        $tot_gs++;
        shift(@p); shift(@p);
        foreach my $g (@p) { $all_genes{$g}++; }

        print OUT "$_\n";
    }
}

my $tot_genes = scalar(keys %all_genes);

print "\nNo. genesets: $tot_gs\nNo. genes: $tot_genes\n\n";

close OUT;
