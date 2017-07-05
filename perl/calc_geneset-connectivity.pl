#!/usr/bin/perl

# ================================================
# Name : calc_geneset-connectivity.pl
# Purpose : Calculate connectivity of network genes to a given geneset
# Created : 01-05-2014
# Last Modified : Wed 25 Mar 2015 10:54:40 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(shuffle);

sub get_cvint;

my $dat2dab = '/Genomics/fgrid/function/sleipnir-build/bin/Dat2Dab';


my ($help, $idab, $igenew, $igenew2, $iseed, $igmt, $ibgg, $igdeg, $iexe);
my $iming = 5; my $imaxg = 200;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(   'help' => \$help,
            'idab=s' => \$idab,
          'igenew=s' => \$igenew,
         'igenew2=s' => \$igenew2,
            'iexe=s' => \$iexe,
           'iseed=s' => \$iseed,
            'igmt=s' => \$igmt,
            'ibgg=s' => \$ibgg,
           'iming=i' => \$iming,
           'imaxg=i' => \$imaxg,
           'igdeg=s' => \$igdeg    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


(my $dtag = $idab) =~ s/^.*\///g; $dtag =~ s/\.da[tb]$//g;
(my $stag = $iseed) =~ s/\.genes$//g; $stag =~ s/^.*\///g;

# Get all edges incident on the seed genes
my $odat = $stag.'.'.$dtag.'.dat';
unless(-e $odat) {
    `$dat2dab -i $idab -D $iseed -o $odat`; }


# In case background-correction is needed, index gene degree
my %gdeg = (); my $ngenes;
if($igdeg) {
    open DEG, "$igdeg" or die "Can't open $igdeg!";
    while (<DEG>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $gdeg{$p[0]} = $p[1]; }
    close DEG;

    $ngenes = scalar keys %gdeg; }


# Index seed genes
my %seed = ();
open SET, "$iseed" or die "Can't open $iseed!";
while (<SET>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $seed{$_}++; }
close SET;

my @seeda = shuffle keys %seed;
my $nseed = scalar @seeda;


# Index target genes in the network (those that are not seeds)
unless($igenew) {
    ($igenew = $idab) =~ s/\.da[tb]$/\.genes/g; $igenew =~ s/^.*\///g;
    unless(-e $igenew) {
        `$dat2dab -i $idab -E > $igenew`; } }

my %target = my %gene_wgt = ();
open GEW, "$igenew" or die "Can't open $igenew!";
while (<GEW>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $gene_wgt{$p[0]} = 1;
    if($#p > 0) { $gene_wgt{$p[0]} = $p[1]; }
    if(exists $seed{$p[0]}) { next; }
    $target{$p[0]}++; }
close GEW;

# Randomly split target genes into groups, each group corresponding to one seed
# gene. This allows calculation of a final gene score based on the same number
# of seed genes, irrespective of wheather the gene is a source or a target
my %cvint = get_cvint(\%target, $nseed);
my %seed_target = ();
for(my $i=0; $i<=$#seeda; $i++) {
    foreach my $g (keys %{$cvint{$i}}) {
        $seed_target{$seeda[$i]}{$g}++; } }


# Index edges to be excluded, if needed
my %ex_edges = ();
if($iexe) {
    open EXE, "$iexe" or die "Can't open $iexe!";
    while (<EXE>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        my $e = join '__', sort($p[0], $p[1]);
        $ex_edges{$e}++; }
    close EXE; }


# Calcualte gene connectivity score
#my %gconn = my %count = ();
my %numr = my %deno = ();
open DAT, "$odat" or die "Can't open $odat!";
while (<DAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $e = join '__', sort($p[0], $p[1]);
    if($iexe and (exists $ex_edges{$e})) { next; }

    my $edge_wgt = $p[2];

    if(exists $seed{$p[0]}) {
        my $a = $p[0]; my $b = $p[1];
        if(($nseed > 1) and (exists $seed_target{$a}{$b})) { next; }

        #$edge_wgt *= $gene_wgt{$b};
        #if($iwseed) {
        #    $edge_wgt *= $gene_wgt{$a}; }

        $numr{$b} += ($edge_wgt * $gene_wgt{$a} * $gene_wgt{$b});
        $deno{$b} ++; }

    #unless(exists $gconn{$b}) { $gconn{$b} = 0; }
    #$count{$b}++;
    #$gconn{$b} += (($edge_wgt - $gconn{$b}) / $count{$b}); }
 
    if(exists $seed{$p[1]}) {
        my $a = $p[1]; my $b = $p[0];
        if(($nseed > 1) and (exists $seed_target{$a}{$b})) { next; }

        #$edge_wgt *= $gene_wgt{$b};
        #if($iwseed) {
        #    $edge_wgt *= $gene_wgt{$a}; }

        $numr{$b} += ($edge_wgt * $gene_wgt{$a} * $gene_wgt{$b});
        $deno{$b} ++; } }

#unless(exists $gconn{$b}) { $gconn{$b} = 0; }
#$count{$b}++;
#$gconn{$b} += (($edge_wgt - $gconn{$b}) / $count{$b}); } }
close DAT;


if($nseed == 1) {
    my @rg = shuffle keys %numr;
    $numr{$seeda[0]} = $numr{$rg[0]};
    $deno{$seeda[0]} = $deno{$rg[0]}; }


# Index background genes for target genesets
my %all_gs_genes = ();
if($ibgg) {
    open BGG, "$ibgg" or die "Can't open $ibgg!";
    while (<BGG>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $all_gs_genes{$_}++; }
    close BGG; }


my %gs_genes = my %gs_desc = my %gs_size = ();
if($igmt) {
    open GMT, "$igmt" or die "Can't open $igmt!";
    while (<GMT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        my $gs = shift @p;
        $p[0] =~ s/ \([0-9]*\)$//g; $gs_desc{$gs} = shift @p;

        foreach my $g (@p) {
            unless(exists $numr{$g}) { next; }
            unless($ibgg) { $all_gs_genes{$g}++; }
            $gs_genes{$gs}{$g}++; }

        $gs_size{$gs} = scalar keys %{$gs_genes{$gs}};

        if(($gs_size{$gs} < $iming) or ($gs_size{$gs} > $imaxg)) {
            delete $gs_genes{$gs}; } } }


my %gene_score = (); my $mean = my $sd = my $n = 0;
# Print gene connectivity scores
my $otab1 = $stag.'.'.$dtag.'_genes.conn';
open TAB1, ">$otab1";
foreach my $g (sort keys %numr) {
    my $v = ( $numr{$g} / $deno{$g} );
    if($igdeg) {
        $v /= ( $gdeg{$g} / $ngenes ); }

    if($igmt) {
        unless(exists $all_gs_genes{$g}) { next; }
        $gene_score{$g} = $v; }

    $n++;
    my $om = $mean;
    $mean += ($v - $om) / $n;
    $sd += ($v - $om)*($v - $mean);

    print TAB1 "$g\t", sprintf("%.6f\n", $v); }

print TAB1 "#mean\t", sprintf("%.6f\n", $mean);
print TAB1 "#sd\t", sprintf("%.6f\n", sqrt($sd/$n));


if($igmt) {
    (my $otab2 = $otab1) =~ s/_genes/_gs/g;
    open TAB2, ">$otab2";
    print TAB2 "#id\tdesc\tsize\tscore\n";

    foreach my $gs (keys %gs_genes) {
        my $gs_score = 0;
        foreach my $g (keys %{$gs_genes{$gs}}) {
            $gs_score += $gene_score{$g}; }

        $gs_score /= $gs_size{$gs};
        print TAB2 "$gs\t$gs_desc{$gs}\t$gs_size{$gs}\t";
        printf TAB2 "%.6f\n", $gs_score; }
    close TAB2; }
close TAB1;





# Given instances and size of cv-fold, assigns instances to cv-fold
sub get_cvint {
    my $gen = shift;
    my $cvk = shift;

    # Randomize the instances
    my @rnd = shuffle keys %$gen;

    my %cvi = (); my $i = 0;

    foreach (@rnd) {
        $cvi{$i++ % $cvk}{$_}++; }

    return %cvi; }

