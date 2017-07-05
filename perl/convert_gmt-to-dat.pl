#!/usr/bin/perl

# ================================================
# Name : conver_gmt-to-dat.pl
# Purpose : Convert genesets to gold-std positive edges
# Created : 16-05-2014
# Last Modified : Fri 16 May 2014 10:55:24 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

my ($igmt, $imaxg, $odat) = @ARGV;


my %edges = (); my $ngs = 0;
open GMT, "$igmt" or die "Can't open $igmt!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    if((($#p-1) > $imaxg) or ($#p == 2)) { next; }

    $ngs++;
    shift @p; shift @p;
    for(my $i=0; $i<$#p; $i++) {
        for(my $j=($i+1); $j<=$#p; $j++) {
            my $e = join '__', sort($p[$i], $p[$j]);
            $edges{$e}++; } } }
close GMT;

print "\n$ngs genesets\n";


my $ne = 0;
open DAT, ">$odat";
foreach my $e (sort keys %edges) {
    $ne++;
    my ($g1, $g2) = split '__', $e;
    print DAT "$g1\t$g2\t1\n"; }
close DAT;

print "$ne edges\n\n";
