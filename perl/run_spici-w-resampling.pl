#!/usr/bin/perl

# ================================================
# Name : run_spici-w-resampling.pl
# Purpose : Run spici clustering with edge-resampling
# Created : 24-12-2014
# Last Modified : Tue 06 Jan 2015 12:00:04 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Graph::Undirected;

my $spici = '/home/arjunk/software/clustering/SPICi/bin/spici';


my $idat = $ARGV[0];
my $iden = $ARGV[1];
my $odat = $ARGV[2];
my $ogmt = $ARGV[3];

my $nitr = 100;


chomp(my $l = `wc -l $idat`); $l =~ s/ .*$//g;
my $sl = int(4*$l/5);
#my $sl = int(2*$l/3);

#my %itr_gene = my %itr_edge = my %edges = ();
my %itr_gene = my %edges = ();
for(my $i=0; $i<$nitr; $i++) {
    unless($i % 10) { print "\nitr $i ..."; }

    `perl -MList::Util -e 'print List::Util::shuffle <>' $idat | head -n $sl > subsamp.dat`;
    #`Dat2Dab -i $idat -u 0.80 -o subsamp.dat`;
    `Dat2Dab -i subsamp.dat -E > subsamp.genes`;
    `$spici -i subsamp.dat -d $iden -m 2 -o subsamp.clu`;

    open GEN, "subsamp.genes" or die "Can't open subsamp.genes!";
    while (<GEN>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $itr_gene{$i}{$_}++; }
    close GEN;

    open CLU, "subsamp.clu" or die "Can't open subsamp.clu!";
    while (<CLU>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;

        for(my $m=0; $m<$#p; $m++) {
            for(my $n=($m+1); $n<=$#p; $n++) {
                my $e = join '__', sort($p[$m], $p[$n]);
                #$itr_edge{$i}{$e}++;
                $edges{$e}++; } } }
    close CLU; }
print "\n"; `rm -f subsamp.*`;


print "\nPrinting DAT ...";
open DAT, ">$odat";

my $icut = 0.8;
my $graph = Graph::Undirected->new;
my $nedges = 0;

foreach my $e (keys %edges) {
    $nedges++;
    my ($g1, $g2) = split '__', $e;

    my $nincl = 0;
    for(my $i=0; $i<$nitr; $i++) {
        if((exists $itr_gene{$i}{$g1}) and (exists $itr_gene{$i}{$g2})) {
            $nincl++; } }

    if($nincl == 0) { die "\n$nedges\t$g1\t$g2\t$edges{$e}\t$nincl\n"; }
    my $score = $edges{$e}/$nincl;
    print DAT "$g1\t$g2\t", sprintf("%.6f\n", $score);

    if($score < $icut) { next; }
    $graph->add_edge($g1, $g2); }
close DAT;


print "\nPrinting GMT ...";
open GMT, ">$ogmt";
my @cc = $graph->connected_components();
my $idx = 0;
foreach my $c (@cc) {
    $idx++;
    my $cid = 'C'.sprintf("%04d", $idx);

    print GMT "$cid\t$cid (", scalar @$c;
    print GMT ")\t", join "\t", @$c;
    print GMT "\n"; }
close GMT;

print "\nDONE\n\n";

