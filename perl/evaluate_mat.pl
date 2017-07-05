#!/usr/bin/perl

# ================================================
# Name : eval_mat.pl
# Purpose : Evaluate genelist again matrix of predictions
# Created : 02-06-2014
# Last Modified : Mon 02 Jun 2014 01:40:44 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

sub evaluate;

my ($igenes, $ibg, $ipred, $oeval) = @ARGV;


my %gene_truel = ();
open BG, "$ibg" or die "Can't open $ibg!";
while (<BG>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $gene_truel{$_} = -1; }
close BG;


open GP, "$igenes" or die "Can't open $igenes!";
while (<GP>) {
    if($_ =~ /^#/) { next; }
    chomp($_); unless(exists $gene_truel{$_}) { next; }
    $gene_truel{$_} = 1; }
close GP;


my %gene_preds = (); my $l = 0; my @meth = ();
open PR, "$ipred" or die "Can't open $ipred!";
while (<PR>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    if($l == 0) {
        @meth = @p; $l++; next; }

    unless($gene_truel{$p[0]}) { next; }

    for(my $j=1; $j<=$#p; $j++) {
        if($p[$j] =~ /NA/) { next; }
        $gene_preds{$meth[$j]}{$p[0]} = $p[$j]; } }
close PR;


open my $EVAL, ">$oeval";
print $EVAL "method\tp\tn\tprior\tauc\tp10r\tp20r\tp50r\tauprc\twauprc\n";
foreach my $m (sort keys %gene_preds) {
    print "$m\n"; print $EVAL "$m";
    evaluate(\%gene_truel, $gene_preds{$m}, $EVAL); }
close $EVAL;


# Evaluate
sub evaluate {
    my $truel = shift;
    my $preds = shift;
    my $FH = shift;

    my $P = my $N = 0;
    foreach (keys %$truel) {
        unless(exists $preds->{$_}) { next; }
        if($truel->{$_} == 1) { $P++; }
        elsif($truel->{$_} == -1) { $N++; } }

    my $auc = my $auprc = my $wauprc = my $p10r = my $p20r = my $p50r = 0;
    my $tp = my $fp = my $tn = my $fn = 0;
    my $in10 = my $in20 = my $in50 = 1;
    my $pr = my $rc = my $w = my $sumw = 0;
    
    foreach (sort {$preds->{$b} <=> $preds->{$a}} keys %{$preds}) {
        unless(exists $truel->{$_}) { next; }

        if($truel->{$_} == 1) {
            $tp++;
            $pr = $tp / ($tp + $fp);
            $rc = $tp / $P;

            if(($rc >= 0.10) and $in10) {
                $p10r = $pr; $in10 = 0; }
            if(($rc >= 0.20) and $in20) {
                $p20r = $pr; $in20 = 0; }
            if(($rc >= 0.50) and $in50) {
                $p50r = $pr; $in50 = 0; }

            $auprc += $pr;

            # $w = $p**($tp-1);
            $w = (-1)*log($rc) / log(2);
            $sumw += $w;
            $wauprc += $w*$pr; }

        elsif($truel->{$_} == -1) {
            $fp++;
            $pr = $tp / ($tp + $fp);
            $auc += $tp; }

        $tn = $N - $fp; $fn = $P - $tp; }

    if(($tp == 0) or ($fp == 0)) {
        print STDERR "warning: Too few +ve true labels or -ve true labels\n";
        $auc = $auprc = $wauprc = 0; }
    else {
        $auc = $auc / $tp / $fp;
        $auprc = $auprc / $tp;
        $wauprc = $wauprc / $sumw; }

#     print $FH "#P\t$P\n#N\t$N\n";
#     print $FH "#PRIOR\t", sprintf("%.3f\n", ($P/($P+$N)));
#     print $FH "#AUC\t", sprintf("%.3f\n", $auc);
#     print $FH "#P10R\t", sprintf("%.3f\n", $p10r);
#     print $FH "#P20R\t", sprintf("%.3f\n", $p20r);
#     print $FH "#P50R\t", sprintf("%.3f\n", $p50r);
#     print $FH "#AUPRC\t", sprintf("%.3f\n", $auprc);
#     print $FH "#WAUPRC\t", sprintf("%.3f\n", $wauprc);
    print $FH "\t$P\t$N";
    print $FH "\t", sprintf("%.3f", ($P/($P+$N)));
    print $FH "\t", sprintf("%.3f", $auc);
    print $FH "\t", sprintf("%.3f", $p10r);
    print $FH "\t", sprintf("%.3f", $p20r);
    print $FH "\t", sprintf("%.3f", $p50r);
    print $FH "\t", sprintf("%.3f", $auprc);
    print $FH "\t", sprintf("%.3f\n", $wauprc);
}
