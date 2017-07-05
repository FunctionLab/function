#!/usr/bin/perl

# ================================================
# Name : run_spici.pl
# Purpose : Run spici clustering over a number of density parameters
# Created : 18-12-2014
# Last Modified : Tue 23 Dec 2014 05:25:19 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

use lib '/Genomics/Users/arjunk/software/lib/perl5/';
use PDF;

my $spici = '/home/arjunk/software/clustering/SPICi/bin/spici';


my $idat = $ARGV[0];
my $istd = $ARGV[1];
my $otab = $ARGV[2];


(my $itag = $idat) =~ s/\.dat$//g; $itag =~ s/^.*\///g;


my $igen = $itag.'.genes';
unless(-e $igen) {
    `Dat2Dab -i $idat -E > $igen`; }

my %genes = ();
open GEN, "$igen" or die "Can't open $igen!";
while (<GEN>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $genes{$_}++; }
close GEN;


(my $istdsub = $istd) =~ s/\.da[tb]$/\.sub\.dat/g; $istdsub =~ s/^.*\///g;
`Dat2Dab -i $istd -g $igen -o $istdsub`;


my %pose = my %nege = my %stdg = ();
open STD, "$istdsub" or die "Can't open $istdsub!";
while (<STD>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    $stdg{$p[0]}++; $stdg{$p[1]}++;
    my $e = join '__', sort($p[0], $p[1]);

    if($p[2] == 1) { $pose{$e}++; }
    elsif($p[2] == 0) { $nege{$e}++; } }
close STD;

my $P = scalar keys %pose;
my $N = scalar keys %nege;


open TAB, ">$otab";
print TAB "#density\tfgenes\tnclu\tnclu.ge5\tbac\tnlogp\n";

#my $m = 0; # sparse graph
#my $m = 1; # dense graph
my $m = 2; # large sparse graph
for(my $d=0.05; $d<1; $d+=0.05) {
    $d = sprintf("%.2f", $d);
    print "$d ...\n";
    my $oclu = $itag.'.d'.$d.'.clust';
    `$spici -i $idat -d $d -m $m -o $oclu`;

    my %clue = (); my $ngen = my $nclu = my $nclu_ge5 = 0;
    open CLU, "$oclu" or die "Can't open $oclu!";
    while (<CLU>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;

        $nclu++; if((scalar @p) >= 5) { $nclu_ge5++; $ngen += scalar @p; }
        for(my $i=0; $i<$#p; $i++) {
            unless(exists $stdg{$p[$i]}) { next; }

            for(my $j=($i+1); $j<=$#p; $j++) {
                unless(exists $stdg{$p[$j]}) { next; }

                my $e = join '__', sort($p[$i], $p[$j]);
                $clue{$e}++; } } }
    close CLU;
    `rm -f $oclu`;

    my $tp = my $tn = 0;
    foreach my $e (keys %pose) {
        if(exists $clue{$e}) { $tp++; } }
    foreach my $e (keys %nege) {
        unless(exists $clue{$e}) { $tn++; } }

    my $bac = 0.5*(($tp/$P) + ($tn/$N));

    #my $lor = log(($tp/(scalar keys %clue))/($P/($P+$N)))/log(2);
    my $pval = hypergeometric_tail(($P+$N), $P, (scalar keys %clue), $tp);

    print TAB "$d\t", sprintf("%.6g", ($ngen/(scalar keys %genes)));
    print TAB "\t$nclu\t$nclu_ge5\t", sprintf("%.6g\t%.6g\n", $bac, (-log($pval)/log(10))); }

print "\n";
