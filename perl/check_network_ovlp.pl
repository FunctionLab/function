#!/usr/bin/perl
use strict;
use warnings;
use PDF;

my($f1, $f2) = @ARGV;
# f1: Network 1 (DAT)
# f2: Network 2 (DAT)
my $par_miss = 0;

open FH, "$f1" or die "Can't open $f1!"; chomp(my @f=<FH>); close FH;
open GH, "$f2" or die "Can't open $f2!"; chomp(my @g=<GH>); close GH;

my (@p, $edge);
my %union_genes = (); my %union_edges = ();
my %common_genes = (); my %common_edges = ();
my %union_posedges = ();

my %net1_genes = (); my %net1_edges = ();
my ($num_net1_posedges, $num_net1_negedges) = (0, 0);
# my %net1_posedges = (); my %net1_negedges = ();
foreach my $e (@f) {
    @p=(); @p=split("\t",$e);
    $edge = join '__', sort($p[0], $p[1]);

    $net1_genes{$p[0]}++; $net1_genes{$p[1]}++;
    $union_genes{$p[0]}++; $union_genes{$p[1]}++;

    $union_edges{$edge}++;
    if(abs($p[2]) > 0) { $net1_edges{$edge}=1; $num_net1_posedges++;
        $union_posedges{$edge}++; }
    else { $net1_edges{$edge}=0; $num_net1_negedges++; }
}

my %net2_genes = (); my %net2_edges = ();
my ($num_net2_posedges, $num_net2_negedges) = (0, 0);
# my %net2_posedges = (); my %net2_negedges = ();
foreach my $e (@g) {
    @p=(); @p=split("\t",$e);
    $edge = join '__', sort($p[0], $p[1]);

    $net2_genes{$p[0]}++; $net2_genes{$p[1]}++;
    $union_genes{$p[0]}++; $union_genes{$p[1]}++;
    if(exists $net1_genes{$p[0]}) { $common_genes{$p[0]}++; }
    if(exists $net1_genes{$p[1]}) { $common_genes{$p[1]}++; }

    $union_edges{$edge}++;
    if(exists $net1_edges{$edge}) { $common_edges{$edge}++; }
    
    if(abs($p[2]) > 0) { $net2_edges{$edge}=1; $num_net2_posedges++;
        $union_posedges{$edge}++; }
    else { $net2_edges{$edge}=0; $num_net2_negedges++; }
}

my $num_net1_genes = scalar(keys %net1_genes);
my $num_net2_genes = scalar(keys %net2_genes);

my $num_common_genes = scalar(keys %common_genes);
my $num_common_edges = scalar(keys %common_edges);
my $num_union_genes = scalar(keys %union_genes);
my $num_union_edges = scalar(keys %union_edges);
my $num_union_posedges = scalar(keys %union_posedges);

my ($num_net1_cposedges, $num_net1_cnegedges) = (0, 0);
my ($num_net2_cposedges, $num_net2_cnegedges) = (0, 0);
my ($num_00_edges, $num_01_edges, $num_10_edges, $num_11_edges) = (0, 0, 0, 0);

foreach my $e (keys %common_edges) {
    if($net1_edges{$e} == 0) {
        $num_net1_cnegedges++;
        if($net2_edges{$e} == 0) { $num_00_edges++; $num_net2_cnegedges++; }
        else { $num_01_edges++; $num_net2_cposedges++; }
    }
    else {
        $num_net1_cposedges++; 
        if($net2_edges{$e} == 0) { $num_10_edges++; $num_net2_cnegedges++; }
        else { $num_11_edges++; $num_net2_cposedges++; }
    }
}

print "\nNet1: $f1\n\tNo. Genes: $num_net1_genes\n";
print "\tNo. pos. edges: $num_net1_posedges\n\tNo. neg. edges: $num_net1_negedges\n";
print "\nNet2: $f2\n\tNo. Genes: $num_net2_genes\n";
print "\tNo. pos. edges: $num_net2_posedges\n\tNo. neg. edges: $num_net2_negedges\n";

print "\nNo. union genes: $num_union_genes\nNo. union edges: $num_union_edges\n";
print "\nNo. common genes: $num_common_genes\nNo. common edges: $num_common_edges\n";

my $prob = hypergeometric_tail($num_union_genes, $num_net1_genes,
    $num_net2_genes, $num_common_genes);
print "Hypergeometric probability: $prob\n";

print "\n\tNet1\t0\t1\nNet2\t$num_common_edges\t$num_net2_cnegedges\t$num_net2_cposedges\n";
print "0\t$num_net1_cnegedges\t$num_00_edges\t$num_01_edges\n";
print "1\t$num_net1_cposedges\t$num_10_edges\t$num_11_edges\n";

$prob = hypergeometric_tail($num_common_edges, $num_net1_cposedges,
    $num_net2_cposedges, $num_11_edges);
print "Hypergeometric probability: $prob\n\n";

