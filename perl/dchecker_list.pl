#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my($ilist, $igmt, $iposfilter, $inegfilter, $otab) = @ARGV;
my (@p, $time, @gsref);

my %gene_score = ();
open FL, "$ilist" or die "Can't open $ilist!";
while(<FL>) {
    chomp($_); @p = split '\t', $_;
    $gene_score{$p[0]} = $p[1];
}
close FL;

@gsref = gs2genes($igmt, \%gene_score);
my %gs_pos = %{$gsref[0]};
my %gs_desc = %{$gsref[1]};
my %gs_size = %{$gsref[2]};
my %all_genes = %{$gsref[3]};

my %gs_subset = ();
open PO, "$iposfilter" or die "Can't open $iposfilter!";
while(<PO>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    $gs_subset{$p[0]}++;
}
close PO;

my %gs_anc = ();
open AC, "$inegfilter" or die "Can't open $inegfilter!";
while(<AC>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    if((exists $gs_subset{$p[0]}) and (exists $gs_subset{$p[1]})) { next; }
    $gs_anc{$p[0]}{$p[1]}++;
}
close AC;

my %gs_neg = ();
foreach my $gs (keys %gs_pos) {
    $gs_neg{$gs} = \%all_genes;
}

sub gs2genes {
    my $file = shift;
    open FH, "$file" or die "Can't open $file!";
        chomp(my @gset=<FH>); close FH;

    my $bgref = shift;

    my (@fl, %genes, %desc, %size, %justgenes, $gs);
    foreach (@gset) {
        @fl = split '\t', $_;
        $gs = shift(@fl);
        $fl[0] =~ s/ \([0-9]*\)$//g;
        $desc{$gs} = shift(@fl);

        foreach my $g (@fl) {
            unless(exists $bgref->{$g}) { next; }
            $genes{$gs}{$g}++;
            $justgenes{$g}++;
        }

        $size{$gs} = scalar keys %{$genes{$gs}};
    }

    return (\%gs_genes, \%gs_desc, \%gs_size, \%justgenes);
}

sub dcheck {
    my $predl = shift; # Original_Class -> Prediction_Score
    my $predt = shift; # File to write out dcheker table
    my $FH = shift; # File handle to write out evaluation summary

    open PL, "$predl"; chomp(my @pred=<PL>); close PL;
    open PT, ">$predt";

    my $tp = my $fp = my $tn = my $fn = 0;
    my $P = my $N = my $pr = my $rc = 0;
    my $auc = my $auprc = my $p10r = my $p20r = my $p50r = 0;
    my $idx = 0; my %pred_score = (); my %true_label = ();
    my $in10 = my $in20 = my $in50 = 1;

    print PT "#Label\tScore\tTP\tFP\tTN\tFN\tPR\tRC\n";

    foreach (@pred) {
        @p = split '\t', $_;
        if($p[0] == 1) { $P++; }
        else { $N++; }
        $true_label{$idx} = $p[0];
        $pred_score{$idx} = $p[1];
        $idx++;
    }

    foreach (sort {$pred_score{$b} <=> $pred_score{$a}} keys %pred_score) {
        if($true_label{$_} == 1) {
            $tp++; $rc = $tp / $P;

            if(($rc >= 0.10) and $in10) {
                $p10r = $tp / ($tp + $fp); $in10 = 0; }
            if(($rc >= 0.20) and $in20) {
                $p20r = $tp / ($tp + $fp); $in20 = 0; }
            if(($rc >= 0.50) and $in50) {
                $p50r = $tp / ($tp + $fp); $in50 = 0; }

            $auprc += $tp / ($tp + $fp); }
        else {
            $fp++;
            $auc += $tp; }

        $tn = $N - $fp; $fn = $P - $tp;
        $pr = $tp / ($tp + $fp); $rc = $tp / ($tp + $fn);

        print PT "$true_label{$_}\t$pred_score{$_}\t$tp\t$fp\t$tn\t$fn\t$pr\t$rc\n";
    }

    if(($tp == 0) or ($fp == 0)) {
        print "warning: Too few +ve true labels or -ve true labels\n";
        $auc = 0; $auprc = 0;
    }
    else {
        $auc = $auc / $tp / $fp;
        $auprc = $auprc / $tp;
    }

    $predl =~ s/\.pred\.labels//g;
    print $FH "$predl\t$P\t$N\t", sprintf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
        $auc, $p10r, $p20r, $p50r, $auprc), "\n";

    close PL; close PT;

    `rm -f $pred_labels`;
}

