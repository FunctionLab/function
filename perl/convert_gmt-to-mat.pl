#!/usr/bin/perl

# ================================================
# Name : convert_gmt-to-mat.pl
# Purpose : Convert GMT to MAT format
# Created : 06-05-2014
# Last Modified : Tue 06 May 2014 12:33:27 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($igmt, $omat) = @ARGV;

open GMT, "$igmt" or die "Can't open $igmt!";
open MAT, ">$omat";

my %gene_gs = my @gs = ();
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $s = shift @p; shift @p;
    push(@gs, $s); print MAT "\t$s";

    foreach my $g (@p) {
        $gene_gs{$g}{$s}++; } }
close GMT;

print MAT "\n";
foreach my $g (sort keys %gene_gs) {
    print MAT "$g";
    foreach my $s (@gs) {
        my $v = 0; if(exists $gene_gs{$g}{$s}) { $v = 1; }
        print MAT "\t$v"; }
    print MAT "\n"; }

