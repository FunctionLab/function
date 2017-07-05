#!/usr/bin/perl

# ================================================
# Name : extract_edges-relevant-to-geneset.pl
# Purpose : Extract links relevant to a particular geneset
# Created : 18-05-2014
# Last Modified : Wed 31 Dec 2014 12:43:01 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

my ($ilink, $iovlp, $igset, $icut, $otag) = @ARGV;

my $osif = $otag.'.sif'; open SIF, ">$osif";
my $onoa = $otag.'.noa.txt'; open NOA, ">$onoa";
my $oleda = $otag.'.link.attrs'; open LEDA, ">$oleda";
my $ooeda = $otag.'.ovlp.attrs'; open OEDA, ">$ooeda";


my %gs_size = my %gs_desc = my %gs_ovlps = ();
open OVLP, "$iovlp" or die "Can't open $iovlp!";
while (<OVLP>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    $gs_size{$p[0]} = $p[1]; $gs_desc{$p[0]} = $p[2];
    $gs_size{$p[3]} = $p[4]; $gs_desc{$p[3]} = $p[5];

    $gs_ovlps{$p[0]}{$p[3]} = $p[11];
    $gs_ovlps{$p[3]}{$p[0]} = $p[11]; }
close OVLP;


my %qry_neighbors = my %gs_links = ();

open LINK, "$ilink" or die "Can't open $ilink!";
while (<LINK>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $u = $p[0]; my $v = $p[3];
    # my $l = $p[12]; my $o = $gs_ovlps{$u}{$v};

    if($p[12] < $icut) { next; }

    $gs_links{$u}{$v} = $p[12];
    $gs_links{$v}{$u} = $p[12];

    if(($u eq $igset) or ($v eq $igset)) {
        $qry_neighbors{$u}++;
        $qry_neighbors{$v}++; } }
close LINK;

#         if($p[12] >= $icut) {
#             print SIF "$p[0]\tlink\t$p[3]\n";
#             print LEDA "$p[0] (link) $p[3] = $p[12]\n"; }
# 
#         print "$p[0]\t$p[3]\t$p[8]\n";
#         if($p[8] >= $icut) {
#             print SIF "$p[0]\tovlp\t$p[3]\n";
#             print OEDA "$p[0] (ovlp) $p[3] = $p[8]\n"; }
# 
#         if($p[0] eq $igset) {
#             $qry_neighbors{$p[3]}++; }
#         
#         elsif($p[3] eq $igset) {
#             $qry_neighbors{$p[0]}++; } } }


print NOA "node\tsize\tname\n";
# print NOA "$igset\t$gs_size{$igset}\t$gs_desc{$igset}\n";

my @neigh = sort keys %qry_neighbors;
foreach my $ngs (@neigh) {
    print NOA "$ngs\t$gs_size{$ngs}\t$gs_desc{$ngs}\n"; }


print "\n";
for(my $i=0; $i<$#neigh; $i++) {
    my $u = $neigh[$i];
    for(my $j=($i+1); $j<=$#neigh; $j++) {
        my $v = $neigh[$j];

        if(exists $gs_ovlps{$u}{$v}) {
            my $o = $gs_ovlps{$u}{$v};
            if($o >= $icut) {
                print SIF "$u\tovlp\t$v\n";
                print OEDA "$u (ovlp) $v = $o\n"; } }
        
        if(exists $gs_links{$u}{$v}) {
            my $l = $gs_links{$u}{$v};
            if($l >= $icut) {
                print SIF "$u\tlink\t$v\n";
                print LEDA "$u (link) $v = $l\n"; } } } }

close SIF; close NOA;
close LEDA; close OEDA;


