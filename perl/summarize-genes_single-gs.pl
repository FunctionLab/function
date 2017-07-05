#!/usr/bin/perl

# ================================================
# Name : summarize_genes-to-genesets.pl
# Purpose : Summarize gene-scores to geneset-scores
# Created : 01-05-2014
# Last Modified : Tue 10 Feb 2015 04:41:40 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(shuffle);


my ($help, $itab, $igs, $iing, $otxt);
my $irand = 25000;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(  'help' => \$help,
           'itab=s' => \$itab,
            'igs=s' => \$igs,
           'iing=s' => \$iing,
          'irand=i' => \$irand,
           'otxt=s' => \$otxt   ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %gene_score = ();
open TAB, "$itab" or die "Can't open $itab!";
while (<TAB>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    if($p[1] eq 'NA') { next; }
    $gene_score{$p[0]} = $p[1]; }
close TAB;


my %all_genes = ();
if($iing) {
    open ING, "$iing" or die "Can't open $iing!";
    while (<ING>) {
        if($_ =~ /^#/) { next; }
        chomp($_);
        unless(exists $gene_score{$_}) { next; }
        $all_genes{$_}++; }
    close ING; }
else {
    map { $all_genes{$_}++ } keys %gene_score; }

my @all_gene_scores = (); my $mean_gene_score = my $num_gene_score = 0;
foreach my $g (keys %all_genes) {
    $num_gene_score++;
    $mean_gene_score += ($gene_score{$g} - $mean_gene_score) / $num_gene_score;
    push(@all_gene_scores, $gene_score{$g}); }


my $gs_score = my $gs_size = 0;
open GS, "$igs" or die "Can't open $igs!";
while (<GS>) {
    if($_ =~ /^#/) { next; }
    chomp($_);

    if($iing) {
        unless(exists $all_genes{$_}) { next; } }

    $gs_size++;
    $gs_score += (($gene_score{$_} - $gs_score) / $gs_size); }
close GS;


my $r = my $exc = 1; my $gs_pvalue;
do {
    my @shuf_idx = shuffle(0..$#all_gene_scores);
    my @pick_idx = @shuf_idx[ 0 .. ($gs_size-1) ];
    my @rands = @all_gene_scores[ @pick_idx ];

    my $rmean = 0;
    for(my $s=0; $s<$gs_size; $s++) {
        $rmean += (($rands[$s] - $rmean) / ($s + 1)); }

    if($rmean >= $gs_score) {
        $exc++; }

    $gs_pvalue = ( $exc / ($r+1) );
    $r++;
} until (( $exc > 10 ) or ($r == $irand));


open OUT, ">$otxt";
print OUT "#gs\tsize\tscore\tpvalue\n";
print OUT "$igs\t$gs_size\t";
printf OUT "%.6g\t%.6g\n",$gs_score, $gs_pvalue;
close OUT;


# Get P-value based on empirical distribution

