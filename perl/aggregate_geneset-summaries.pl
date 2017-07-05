#!/usr/bin/perl

# ================================================
# Name : aggregate_geneset-summaries.pl
# Purpose : Aggreagte geneset-summaries into geneset-associations
# Created : 13-05-2014
# Last Modified : Sun 25 Jan 2015 09:18:33 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Statistics::Distributions;

sub get_esfdr;


my ($help, @isumm, $iovlp, $otab);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(   'help' => \$help,
        'isumm=s{,}' => \@isumm,
           'iovlp=s' => \$iovlp,
            'otab=s' => \$otab    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %gs_desc = my %gs_size = ();
my %gs_ovlp_ncomm = my %gs_ovlp_pvalue = my %gs_ovlp_qvalue = ();

open OVLP, "$iovlp" or die "Can't open $iovlp!";
while (<OVLP>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    $gs_size{$p[0]} = $p[1]; $gs_desc{$p[0]} = $p[2];
    $gs_size{$p[3]} = $p[4]; $gs_desc{$p[3]} = $p[5];

    my $tag = join '__', sort($p[0], $p[3]);

    $gs_ovlp_ncomm{$tag} = $p[6];
    $gs_ovlp_pvalue{$tag} = $p[10];
    $gs_ovlp_qvalue{$tag} = $p[11]; }
close OVLP;


my %link_mean = my %link_pvalue = ();

foreach my $summ (@isumm) {
    (my $gs = $summ) =~ s/\.summ$//g;
    $gs =~ s/__/:/g;

    open SUMM, "$summ" or die "Can't open $summ!";
    while (<SUMM>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        if($p[4] == 1) { next; }

        $link_mean{$gs}{$p[0]} = $p[3];
        $link_pvalue{$gs}{$p[0]} = $p[4]; }
    close SUMM; }


my %gs_link_mean = my %gs_link_pvalue = my %gs_link_zscore = ();

foreach my $gs_tag (keys %gs_ovlp_ncomm) {
    my ($gs1, $gs2) = split '__', $gs_tag;
    unless((exists $link_mean{$gs1}) and (exists $link_mean{$gs2})) {
        next; }
    unless((exists $link_mean{$gs1}{$gs2}) and (exists $link_mean{$gs2}{$gs1})) {
        next; }

    $gs_link_mean{$gs_tag} = ( $link_mean{$gs1}{$gs2} + $link_mean{$gs2}{$gs1} ) / 2;

    if(($link_pvalue{$gs1}{$gs2} == 1) or ($link_pvalue{$gs1}{$gs2} == 0)) {
        die "\t$gs1\t$gs2\n"; }
    if(($link_pvalue{$gs2}{$gs1} == 1) or ($link_pvalue{$gs2}{$gs1} == 0)) {
        die "\t$gs2\t$gs1\n"; }

    my $z12 = Statistics::Distributions::udistr($link_pvalue{$gs1}{$gs2});
    my $z21 = Statistics::Distributions::udistr($link_pvalue{$gs2}{$gs1});

    $gs_link_zscore{$gs_tag} = ( $z12 + $z21 ) / sqrt(2);
    $gs_link_pvalue{$gs_tag} = Statistics::Distributions::uprob($gs_link_zscore{$gs_tag}); }


my %gs_link_esq = get_esfdr(\%gs_link_pvalue);
#my %gs_link_bgz = bgcor(\%gs_link_zscore);

open TAB, ">$otab";
# print TAB "#gs1\tgs1.size\tgs1.desc\tgs2\tgs2.size\tgs2.desc\tovlp.size\tovlp.pval\tovlp.esq\tlink.mean\tlink.pval\tlink.zscore\tlink.esq\tlink.bgz\n";
print TAB "#gs1\tgs1.size\tgs1.desc\tgs2\tgs2.size\tgs2.desc\tovlp.size\tovlp.pval\tovlp.esq\tlink.mean\tlink.pval\tlink.zscore\tlink.esq\n";
foreach my $gs_tag (sort {$gs_link_esq{$b} <=> $gs_link_esq{$a}} keys %gs_link_pvalue) {
    my ($gs1, $gs2) = split '__', $gs_tag;

    print TAB "$gs1\t$gs_size{$gs1}\t$gs_desc{$gs1}\t";
    print TAB "$gs2\t$gs_size{$gs2}\t$gs_desc{$gs2}\t";
    print TAB "$gs_ovlp_ncomm{$gs_tag}\t$gs_ovlp_pvalue{$gs_tag}\t$gs_ovlp_qvalue{$gs_tag}\t";
    printf TAB "%.6f\t%.6g\t", $gs_link_mean{$gs_tag}, $gs_link_pvalue{$gs_tag};
    printf TAB "%.6f\t%.6f\t", $gs_link_zscore{$gs_tag}, $gs_link_esq{$gs_tag};
    printf TAB "\n"; }
    #if(exists $gs_link_bgz{$gs_tag}) {
    #    printf TAB "%.6f\n", $gs_link_bgz{$gs_tag}; }
    #else { print TAB "NA\n"; } }
close TAB;


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


# Background-correct links
sub bgcor {
    my $link = shift;

    my %node_count = my %node_mean = my %node_sd = (); my $om;
    foreach my $l (keys %$link) {
        my ($a, $b) = split '__', $l;

        $node_count{$a}++;
        $om = 0; if(exists $node_mean{$a}) { $om = $node_mean{$a}; }
        $node_mean{$a} += ( ( $link->{$l} - $om ) / $node_count{$a} );
        $node_sd{$a} += ( $link->{$l} - $om )*( $link->{$l} - $node_mean{$a} );

        $node_count{$b}++;
        $om = 0; if(exists $node_mean{$b}) { $om = $node_mean{$b}; }
        $node_mean{$b} += ( ( $link->{$l} - $om ) / $node_count{$b} );
        $node_sd{$b} += ( $link->{$l} - $om )*( $link->{$l} - $node_mean{$b} ); }

    foreach my $n (keys %node_sd) {
        if($node_count{$n} <= 1) {
            print "$n\t$node_count{$n}\t$node_mean{$n}\n"; next; }
        $node_sd{$n} = sqrt( $node_sd{$n} / ($node_count{$n}-1) ); }

    my %zscore = ();
    foreach my $l (keys %$link) {
        my ($a, $b) = split '__', $l;
        if(($node_count{$a} <= 1) or ($node_count{$b} <= 1)) { next; }
        if(($node_sd{$a} == 0) or ($node_sd{$b} == 0)) {
            next; }
        # die "$l\t$node_count{$a}\t$node_count{$b}\t$node_sd{$a}\t$node_sd{$b}\n"; }

        my $za = ($link->{$l} - $node_mean{$a}) / $node_sd{$a};
        my $zb = ($link->{$l} - $node_mean{$b}) / $node_sd{$b};
        my $z = ( $za + $zb ) / sqrt(2);

        $zscore{$l} = $z; }

    return %zscore; }
