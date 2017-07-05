#!/usr/bin/perl
use strict;
use warnings;

my($in_net1, $in_net2, $out_tab) = @ARGV;
open FH, "$in_net1" or die "Can't open $in_net1!";
    chomp(my @net1=<FH>); close FH;
open GH, "$in_net2" or die "Can't open $in_net2!";
    chomp(my @net2=<GH>); close GH;
open HH, ">$out_tab";

my %net2_edges = (); my (@p, $e);
foreach (@net2) {
    @p = split '\t', $_;
    $e = join '__', sort($p[0], $p[1]);
    $net2_edges{$e} = $p[2];
}

foreach (@net1) {
    @p = split '\t', $_;
    $e = join '__', sort($p[0], $p[1]);
    unless(exists $net2_edges{$e}) { next; }

    print HH "$p[0]\t$p[1]\t$p[2]\t$net2_edges{$e}\t";
    printf HH "%.6g\n", ($p[2]*$net2_edges{$e});
}

close HH;
