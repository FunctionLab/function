#!/usr/bin/perl

# ================================================
# Name : calc_geneset-significance.pl
# Purpose : Summarize gene-scores to geneset-scores
# Created : 01-05-2014
# Last Modified : Thu 09 Apr 2015 09:22:34 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::Distributions;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);


my ($help, $igene_scores, $igs_scores, $ibgs, $isize, $otab, $ipara);
my $irand = 100000; my $iminz = 0.25;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(  'help' => \$help,
   'igene_scores=s' => \$igene_scores,
     'igs_scores=s' => \$igs_scores,
          'isize=i' => \$isize,
          'irand=i' => \$irand,
          'iminz=f' => \$iminz,
           'ibgs=s' => \$ibgs,
            'ipara' => \$ipara,
           'otab=s' => \$otab   ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if(($ipara) and ($isize < 10)) {
    print "\nwarning: parametric test is not appropriate for gene-sets of size < 10."; }


my $time = runtime(); print "\n$time: tic $igene_scores $isize";


my %gene_score = my @all_gene_scores = (); my ($mean_gene_score, $sd_gene_score);
open GES, "$igene_scores" or die "Can't open $igene_scores!";
while (<GES>) {
    chomp($_); my @p = split '\t', $_;
    if($p[1] eq 'NA') { next; }
    if($_ =~ /^#/) {
        if($p[0] =~ /mean/) { $mean_gene_score = $p[1]; }
        elsif($p[0] =~ /sd/) { $sd_gene_score = $p[1]; }
        next; }
    $gene_score{$p[0]} = $p[1];
    push(@all_gene_scores, $p[1]); }
close GES;


my %gs_desc = my %gs_score = my %gs_zscore = ();
my %remain_gs = my %gs_exc = my %gs_pvalue = (); my $ngs = 0;
open GSS, "$igs_scores" or die "Can't open $igs_scores!";
while (<GSS>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    unless($p[2] == $isize) { next; }

    my $gs = $p[0];
    $gs_desc{$gs} = $p[1];
    $gs_score{$gs} = $p[3];
    $gs_zscore{$gs} = ($p[3] - $mean_gene_score) / ($sd_gene_score/sqrt($isize));
    if($ipara) {
        $gs_pvalue{$gs} = sprintf("%.6g", (Statistics::Distributions::uprob($gs_zscore{$gs})));
        next; }

    if($gs_zscore{$gs} < $iminz) {
        next; }

    $remain_gs{$gs}++; $ngs++;
    $gs_exc{$gs} = 1;
    $gs_pvalue{$gs} = (1/$irand); }
close GSS;
print "\n\ttesting $ngs genesets ...";
if($ngs == 0) { exit; }


unless($ipara) {
    if($ibgs) {
        my $r = 0;
        open BGS, "$ibgs" or die "Can't open $ibgs!";
        while (<BGS>) {
            if($_ =~ /^#/) { next; }
            chomp($_); my @p = split '\t', $_;
            if((scalar @p) != $isize) { die "\nrandom gs not same size!\n\n"; }
            $r++;

            my $rsum = 0; foreach my $g (@p) {
                unless(exists $gene_score{$g}) { print "\nno score $g!\n"; exit; }
                $rsum += $gene_score{$g}; }
            my $rmean = ($rsum / $isize);

            foreach my $gs (keys %remain_gs) {
                if($rmean >= $gs_score{$gs}) {
                    $gs_exc{$gs}++;
                    $gs_pvalue{$gs} = ( $gs_exc{$gs} / ($r+1) ); }

                if($r < 10) { next; }

                if( $gs_exc{$gs} >  10 ) {
                    delete $remain_gs{$gs}; } }

            if(scalar keys %remain_gs == 0) { last; }
            if($r == $irand) { last; } }
        close BGS; }
    else {
        for(my $r=1; $r<$irand; $r++) {
            unless($r % 5000) { $time = runtime(); print "\n\t$time: $r ..."; }
            my @shuf_idx = shuffle(0..$#all_gene_scores);
            my @pick_idx = @shuf_idx[ 0 .. ($isize-1) ];
            my @rands = @all_gene_scores[ @pick_idx ];

            my $rsum = 0;
            map { $rsum += $_ } @rands;
            #map { $rsum += $_ } @all_gene_scores[ @pick_idx ];
            my $rmean = ($rsum / $isize);

            foreach my $gs (keys %remain_gs) {
                if($rmean >= $gs_score{$gs}) {
                    $gs_exc{$gs}++;
                    $gs_pvalue{$gs} = ( $gs_exc{$gs} / ($r+1) ); }

                if($r < 10) { next; }

                if( $gs_exc{$gs} >  10 ) {
                    delete $remain_gs{$gs}; } }

            if(scalar keys %remain_gs == 0) { last; } } } }

#my %gs_qvalue = get_esfdr(\%gs_pvalue);


open OUT, ">$otab";
print OUT "#gs\tdesc\tsize\tmean.score\tzscore";
#print OUT "\tpvalue\tneglog.qvalue\n";
print OUT "\tpvalue\n";
#foreach my $gs (sort {$gs_qvalue{$b} <=> $gs_qvalue{$a}} keys %gs_qvalue) {
foreach my $gs (sort {$gs_pvalue{$a} <=> $gs_pvalue{$b}} keys %gs_pvalue) {
    print OUT "$gs\t$gs_desc{$gs}\t$isize";
    printf OUT "\t%.6f\t%.6f", $gs_score{$gs}, $gs_zscore{$gs};
    #printf OUT "\t%.6g\t%.6f\n", $gs_pvalue{$gs}, $gs_qvalue{$gs}; }
    printf OUT "\t%.6g\n", $gs_pvalue{$gs}; }
close OUT;


$time = runtime(); print "\n$time: toc $igene_scores $isize\n\n";



# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $pref = shift;
    
    my $minp = 1; my $ntests = 0;
    foreach (keys %{$pref}) {
        $ntests++;
        if(($pref->{$_} != 0) 
                and ($pref->{$_} < $minp)) {
            $minp = $pref->{$_}; } }

    my %esfdr = (); my $rank = $ntests;
    my $prev_qval = my $curr_qval = 1;

    foreach (sort {$pref->{$b} <=> $pref->{$a}} keys %{$pref}) {
        my $pvalue = $pref->{$_};

        if($pvalue == 0) {
            $curr_qval = ($minp * $ntests) / ($rank * 10); }
        else {
            $curr_qval = ($pvalue * $ntests) / $rank; }

        my $qvalue = 1;
        if($prev_qval < $curr_qval) {
            $qvalue = $prev_qval; }
        else {
            $qvalue = $curr_qval;
            $prev_qval = $curr_qval; }

        $esfdr{$_} = (-1)*log($qvalue)/log(10);
        if($esfdr{$_} < 0) { $esfdr{$_} = 0; }

        $rank--; }

    return %esfdr; }

