#!/usr/bin/perl

# ================================================
# Name : summarize_genes-to-genesets.pl
# Purpose : Summarize gene-scores to geneset-scores
# Created : 01-05-2014
# Last Modified : Tue 05 Apr 2016 03:16:37 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::Distributions;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);


my ($help, $itab, $igmt, $otab, $ipara, $igenew, $ibgg);
#my $irand = 25000; my $iming = 5; my $imaxg = 100;
my $irand = 100000; my $iming = 5; my $imaxg = 500; my $iminz = -5;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(  'help' => \$help,
           'itab=s' => \$itab,
         'igenew=s' => \$igenew,
           'igmt=s' => \$igmt,
           'ibgg=s' => \$ibgg,
          'iming=i' => \$iming,
          'imaxg=i' => \$imaxg,
          'iminz=f' => \$iminz,
          'irand=i' => \$irand,
            'ipara' => \$ipara,
           'otab=s' => \$otab   ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if(($ipara) and ($iming < 10)) {
    print "\nwarning: parametric test is not appropriate for gene-sets of size < 10."; }


my $time = runtime(); print "\n$time: tic $itab $igmt";


my %gene_score = my %gene_wgt = ();
open TAB, "$itab" or die "Can't open $itab!";
while (<TAB>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $gene_wgt{$p[0]} = 1;
    if($p[1] eq 'NA') { next; }
    $gene_score{$p[0]} = $p[1]; }
close TAB;


if($igenew) {
    open GEN, "$igenew" or die "Can't open $igenew!";
    while (<GEN>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $gene_wgt{$p[0]} = $p[1]; }
    close GEN; }


my %bg_genes = ();
if($ibgg) {
    open BGG, "$ibgg" or die "Can't open $ibgg!";
    while (<BGG>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $bg_genes{$p[0]}++; }
    close BGG; }


my %gs_desc = my %gs_size = ();
my %size_gs = my %gs_score = ();
my %all_gene_sg = ();

open GMT, "$igmt" or die "Can't open $igmt!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    my $gs = shift @p;
    $p[0] =~ s/ \([0-9]*\)$//g; $gs_desc{$gs} = shift @p;

    my %genes = ();
    foreach my $g (@p) {
        unless(exists $gene_score{$g}) { next; }
        if($ibgg) {
            if(exists $bg_genes{$g}) {
                $all_gene_sg{$g}++; } }
        else {
            $all_gene_sg{$g}++; }
        $genes{$g}++; }

    my $size = scalar keys %genes;
    if(($size < $iming) or ($size > $imaxg)) { next; }

    $gs_size{$gs} = $size;
    $size_gs{$size}{$gs}++;

    my $numr = my $deno = 0;
    foreach my $g (keys %genes) {
        $numr += ($gene_score{$g} * $gene_wgt{$g});
        $deno ++; }

    $gs_score{$gs} = ( $numr / $deno ); }

#$gs_score{$gs} = 0; my $n = 0;
#foreach my $g (keys %genes) {
#$n++;
#$gs_score{$gs} += (($gene_score{$g} - $gs_score{$gs}) / $n); } }
close GMT;


my @all_gene_scores = (); my $mean_gene_score = my $sd_gene_score = my $num_gene_score = 0;
foreach my $g (keys %all_gene_sg) {
    $num_gene_score++;
    my $s = ($gene_score{$g} * $gene_wgt{$g});
    my $om = $mean_gene_score;
    $mean_gene_score += ($s - $om) / $num_gene_score;
    $sd_gene_score += ($s - $om)*($s - $mean_gene_score);
    push(@all_gene_scores, $s); }

$sd_gene_score = sqrt($sd_gene_score / $num_gene_score);
#print "mean gene score: $mean_gene_score\n";


my %gs_zscore = my %gs_pvalue = (); my $nsize = 0;
foreach my $size (sort {$a <=> $b} keys %size_gs) {
    my $ngs = 0;
    $nsize++;
    $time = runtime(); print "\n\t$time:\t$nsize\t$size\t", (scalar keys %{$size_gs{$size}});
    foreach my $gs (keys %{$size_gs{$size}}) {
        #if($gs_score{$gs} <= $mean_gene_score) {
        my $z = (($gs_score{$gs} - $mean_gene_score) / ($sd_gene_score/sqrt($size)));
        $gs_zscore{$gs} = $z;

        if($ipara) {
            $gs_pvalue{$gs} = sprintf("%.6g",
                (Statistics::Distributions::uprob($z))); }
        else {
            if($z < $iminz) {
                delete $size_gs{$size}{$gs}; }
            else {
                $ngs++; } } }

    unless($ipara) {
        print "\t$ngs";
        if($ngs == 0) { next; }

        my $pval = emp_dist(\@all_gene_scores, \%gs_score, $size_gs{$size}, $size);

        foreach my $gs (keys %$pval) {
            $gs_pvalue{$gs} = $pval->{$gs}; } } }


my %gs_qvalue = get_esfdr(\%gs_pvalue);
open OUT, ">$otab";
print OUT "#gs\tdesc\tsize\tmean.score\tzscore";
print OUT "\tpvalue\tneglog.qvalue\n";
foreach my $gs (sort {$gs_qvalue{$b} <=> $gs_qvalue{$a}} keys %gs_qvalue) {
    print OUT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    printf OUT "\t%.6f\t%.6f", $gs_score{$gs}, $gs_zscore{$gs};
    printf OUT "\t%.6g\t%.6f\n", $gs_pvalue{$gs}, $gs_qvalue{$gs}; }
close OUT;


$time = runtime(); print "\n$time: toc $itab $igmt\n\n";


# Get P-values based on empirical distribution
sub emp_dist {
    my $scores = shift;
    my $gsscore = shift;
    my $gsref = shift;
    my $size = shift;

    my %remain_gs = %{$gsref};

    my %exc = my %pvalue = ();
    foreach my $gs (keys %remain_gs) {
        $exc{$gs} = 1; $pvalue{$gs} = (1/$irand); }

    for(my $r=1; $r<$irand; $r++) {
        my @shuf_idx = shuffle(0..$#{$scores});
        my @pick_idx = @shuf_idx[ 0 .. ($size-1) ];
        my @rands = @{$scores}[ @pick_idx ];

        my $rsum = 0;
        map { $rsum += $_ } @rands;
        my $rmean = ($rsum / $size);

        #my $rmean = 0;
        #for(my $s=0; $s<$size; $s++) {
        #    $rmean += (($rands[$s] - $rmean) / ($s + 1)); }

        foreach my $gs (keys %remain_gs) {
            if($rmean >= $gsscore->{$gs}) {
                $exc{$gs}++;
                $pvalue{$gs} = ( $exc{$gs} / ($r+1) ); }

            if($r < 10) { next; }

            if( $exc{$gs} >  10 ) {
                delete $remain_gs{$gs}; } }

        if(scalar keys %remain_gs == 0) { last; } }

    return \%pvalue; }


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

