#!/usr/bin/perl

# ================================================
# Name : run_louvain.pl
# Purpose : Run all the clustering steps
# Created : 16-12-2014
# Last Modified : Mon 20 Feb 2017 06:42:53 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Graph::Undirected;

#my $community = '/data/home/arjunk/software/clustering/Louvain/Community_latest/community';
#my $hierarchy = '/data/home/arjunk/software/clustering/Louvain/Community_latest/hierarchy';
my $community = '/Genomics/ogtr04/arjunk/software/Louvain/community';
my $hierarchy = '/Genomics/ogtr04/arjunk/software/Louvain/hierarchy';


#my $ibin = $ARGV[0];
#my $iwgt = $ARGV[1];
#my $igen = $ARGV[2];
my $idat = $ARGV[0];
my $icut = $ARGV[1];
my $iming = $ARGV[2];


`perl /Genomics/ogtr04/arjunk/projects/human-autism-fln/revision/net-clustering/run_cluster.pl $idat`;

(my $ibin = $idat) =~ s/\.dat/\.bin/g;
(my $iwgt = $idat) =~ s/\.dat/\.weights/g;
(my $igen = $idat) =~ s/\.dat/\.genes/g;


(my $itag = $ibin) =~ s/\.bin$//g; $itag =~ s/^.*\///g;
my $ogmt = $itag.'.louvain.gmt';
my $odat = $itag.'.louvain.dat';


my $idx = 0; my %idx_gene = ();
open GEN, "$igen" or die "Can't open $igen!";
while (<GEN>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $idx_gene{$idx} = $_;
    $idx++; }
close GEN;


print "\n# $itag\n\n";
my $nitr = 100; my %edges = ();
for(my $i=0; $i<$nitr; $i++) {
    print "$i ... ";

    `$community $ibin -l -1 -w $iwgt > rand.tree`;

    chomp(my @clus = `$hierarchy -l 1 rand.tree`);
    my %clu_gene = ();
    foreach (@clus) {
        my @p = split ' ', $_;
        $clu_gene{$p[1]}{$idx_gene{$p[0]}}++; }

    #my $nc = 0;
    foreach my $c (keys %clu_gene) {
        my @a = sort keys %{$clu_gene{$c}};

        #$nc++;
        #if($nc <= 5) { print "$nc\t$c\t", scalar @a, "\n"; }

        for(my $m=0; $m<$#a; $m++) {
            for(my $n=($m+1); $n<=$#a; $n++) {
                my $e = join '__', sort($a[$m], $a[$n]);
                $edges{$e}++; } } }
}

`rm -f rand.tree`;


print "\nPrinting DAT ...";
open DAT, ">$odat";

#my $icut = 0.9;
my $graph = Graph::Undirected->new;
my $nedges = 0;

foreach my $e (keys %edges) {
    $nedges++;
    my ($g1, $g2) = split '__', $e;

    my $score = $edges{$e}/$nitr;
    print DAT "$g1\t$g2\t", sprintf("%.6f\n", $score);

    if($score < $icut) { next; }
    $graph->add_edge($g1, $g2); }
close DAT;


print "\nPrinting GMT ...";
open GMT, ">$ogmt";
my @cc = $graph->connected_components();

my %cc_ng = my %cc_genes = (); #my $iming = 0;
foreach my $c (@cc) {
    my $ng = scalar @$c;
    if($ng < $iming) { next; }
    $cc_ng{$c} = $ng;
    push(@{$cc_genes{$c}}, @$c); }

my @selc = sort {$cc_ng{$b} <=> $cc_ng{$a}} keys %cc_ng;

my $cnum = 0;
foreach my $c (@selc) {
    $cnum++;
    my $cid = 'C'.sprintf("%04d", $cnum);

    print GMT "$cid\t$cid (", $cc_ng{$c};
    print GMT ")\t", join "\t", @{$cc_genes{$c}};
    print GMT "\n"; }
close GMT;

print "\nDONE\n\n";

