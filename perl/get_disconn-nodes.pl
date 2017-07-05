#!/usr/bin/perl

# ================================================
# Name : get_disconn-nodes.pl
# Purpose : Get disconnected nodes
# Created : 11-01-2016
# Last Modified : Tue 12 Jan 2016 04:33:04 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

my $idat = $ARGV[0];
my $igen = $ARGV[1];
my $ogen; if((scalar @ARGV) >2) { $ogen = $ARGV[2]; }
else { ($ogen = $idat) =~ s/\.dat/\.disconn-genes\.txt/g; }


my %all_nodes = ();
open IBG, "$igen" or die "Can't open $igen!";
while (<IBG>) {
    if($_ =~ /^#/) { next; }
    chomp($_); if($_ =~ /^$/) { next; }
    $all_nodes{$_}++; }
close IBG;


my %node_pair = my %node_deg = (); my $nedges = 0;
open DAT, "$idat" or die "Can't open $idat!";
while (<DAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my ($g1, $g2, $w) = split '\t', $_;
    $nedges++;

    $node_pair{$g1}{$g2} = $w;
    $node_pair{$g2}{$g1} = $w;

    $node_deg{$g1}++;
    $node_deg{$g2}++; }
close DAT;

my $nnodes = scalar keys %node_deg;
my $nreste = $nedges;


do {
    my @anode = sort {$node_deg{$b} <=> $node_deg{$a}} keys %node_pair;
    my $g = $anode[0]; #print "$g\t$node_deg{$g}\n";
    #print "\t", join " ", sort keys %{$node_pair{$g}}; print "\n";

    if($node_deg{$g} == 0) { last; }
    else {
        foreach my $g2 (keys %{$node_pair{$g}}) {
            #print "\t$g2\t$node_pair{$g2}{$g}\t$node_deg{$g2}\n";
            delete $node_pair{$g2}{$g};
            $node_deg{$g2}--;
            $nreste--; }
        delete $node_pair{$g}; }
} until($nreste == 0);


my $ndcnodes = 0;
open GEN, ">$ogen";
foreach my $g (keys %all_nodes) {
    if(exists $node_deg{$g}) {
        if(exists $node_pair{$g}) {
            $ndcnodes++;
            print GEN "$g\n"; } }
    else {
        print GEN "$g\n"; } }
close GEN;


print "\nAll Nodes: ", scalar keys %all_nodes; print "\n";
print "\nTot. Edges: $nedges\nTot. Nodes: $nnodes\n";
print "\nNo. disconnected nodes: $ndcnodes\n\n";
