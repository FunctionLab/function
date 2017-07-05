#!/usr/bin/perl

# ================================================
# Name : filter_gmt.pl
# Purpose : Filter GMT based by gene-list
# Created : 10-06-2014
# Last Modified : Fri 02 Dec 2016 12:49:28 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

my $iming = 0; my $imaxg = 2000000;
my ($igmt, $igenes, $ogmt) = @ARGV;


my %genes = ();
open GE, "$igenes" or die "Can't open $igenes!";
while (<GE>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $genes{$p[0]}++; }
close GE;


open OG, ">$ogmt";
open IG, "$igmt" or die "Can't open $igmt!";
while (<IG>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $gs = shift @p;
    (my $desc = shift @p) =~ s/ \([0-9]*\)$//g;

    my %tmpg = ();
    foreach my $g (@p) {
        unless(exists $genes{$g}) { next; }
        $tmpg{$g}++; }

    my $ng = scalar keys %tmpg;
    if(($ng < $iming) or ($ng > $imaxg)) { next; }

    print OG "$gs\t$desc ($ng)";
    foreach my $g (sort keys %tmpg) {
        print OG "\t$g"; }
    print OG "\n"; }
close IG;
close OG;
