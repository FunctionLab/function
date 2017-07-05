#!/usr/bin/perl

# ================================================
# Name : compare_cross-species_matrices.pl
# Purpose : Compare a pair of matrices from two different species
# Created : 10-02-2014
# Last Modified : Mon 17 Mar 2014 10:34:15 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

sub index_list;
sub get_gscore;
sub gene2fam;
sub get_fscore;
sub rank;
sub mean_sd;
sub cor;
sub rbo;

my ($help, $imat1, $imat2, $iortho, $ifeat1, $ifeat2, $itopn, $itrues, $omat, $ofam, $otab); my $imes = 'rbo';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
         'imat1=s' => \$imat1,
         'imat2=s' => \$imat2,
          'imes=s' => \$imes,
         'itopn=i' => \$itopn,
        'iortho=s' => \$iortho,
        'ifeat1=s' => \$ifeat1,
        'ifeat2=s' => \$ifeat2,
        'itrues=s' => \$itrues,
          'omat=s' => \$omat,
          'ofam=s' => \$ofam,
          'otab=s' => \$otab    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my $time = runtime(); print "\n$time: indexing gene-fam and gene-scores ...";
# Indexing gene-family membership and gene-feature scores
my ($gene_fam, $fam_gene) = gene2fam($iortho);

my ($gscore1, $genes1) = get_gscore($imat1, 'sp1');
my ($gscore2, $genes2) = get_gscore($imat2, 'sp2');
print "\n\tsp1: ", scalar keys %$genes1, " genes; ", scalar keys %$gscore1, " features";
print "\n\tsp2: ", scalar keys %$genes2, " genes; ", scalar keys %$gscore2, " features";

my %selc_fam = ();
foreach my $f (keys %$fam_gene) {
    my $n1 = my $n2 = 0;
    foreach my $g (keys %{$fam_gene->{$f}}) {
        if(exists $genes1->{$g}) { $n1++; }
        elsif(exists $genes2->{$g}) { $n2++; } }

    if(($n1 > 0) and ($n2 > 0)) { $selc_fam{$f}++; } }

print "\n\t#selc families: ", scalar keys %selc_fam, "\n";


$time = runtime(); print "\n$time: indexing feat desc and initizalizing mean/sd ...";
# Indexing feature descriptions and initializing mean/sd
my %feat1 = index_list($ifeat1, 'sp1');
my %feat2 = index_list($ifeat2, 'sp2');

my @feat1a = sort keys %$gscore1;
my @feat2a = sort keys %$gscore2;
my %feat_count = my %feat_mean = my %feat_sd = ();

my %fscore1 = my %fgene1 = my %frank1 = my %fmean1 = my %fsd1 = ();
foreach (@feat1a) {
    $feat_count{$_} = $feat_mean{$_} = $feat_sd{$_} = 0;

    ($fscore1{$_}, $fgene1{$_}) = get_fscore($gscore1->{$_});
    if($imes eq 'cor') {
        $frank1{$_} = rank($fscore1{$_});
        ($fmean1{$_}, $fsd1{$_}) = mean_sd($frank1{$_}); } }

open OMAT, ">$omat";
foreach (@feat2a) {
    print OMAT "\t$_";
    $feat_count{$_} = $feat_mean{$_} = $feat_sd{$_} = 0; }
print OMAT "\n";

my %fscore2 = my %fgene2 = my %frank2 = my %fmean2 = my %fsd2 = ();
foreach (@feat2a) {
    print OMAT "\t$feat2{$_}";

    ($fscore2{$_}, $fgene2{$_}) = get_fscore($gscore2->{$_});
    if($imes eq 'cor') {
        $frank2{$_} = rank($fscore2{$_});
        ($fmean2{$_}, $fsd2{$_}) = mean_sd($frank2{$_}); } }
print OMAT "\n";


$time = runtime(); print "\n$time: recording true similarity scores ...";
# Recording true similarity scores
my %feat_trues = ();
open TS, "$itrues" or die "Can't open $itrues!";
while (<TS>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $e = join '__', sort($p[0], $p[1]);
    $feat_trues{$e} = $p[2]; }
close TS;


$time = runtime(); print "\n$time: calculating sim b/w feature-pairs ...";
# Calcualting similarity between all pairs of features
my %featpair_sim = (); my $npair = 0;
foreach my $ft1 (@feat1a) {
    foreach my $ft2 (@feat2a) {
        $npair++;
        unless($npair % 500) {
            $time = runtime(); print "\n\t$time: $npair pairs ..."; }

        my $sim;
        if($imes eq 'rbo') {
            # 0.985, 0.990, 0.995
            $sim = rbo($fscore1{$ft1}, $fscore2{$ft2}, 0.99); }
        elsif($imes eq 'cor') {
            # $sim = cor($fscore1{$ft1}, $fscore2{$ft2}); }
            $sim = cor($frank1{$ft1}, $fmean1{$ft1}, $fsd1{$ft1}, $frank2{$ft2}, $fmean2{$ft2}, $fsd2{$ft2}); }

        $feat_count{$ft1}++;
        my $om = $feat_mean{$ft1};
        $feat_mean{$ft1} += (($sim - $om) / $feat_count{$ft1});
        $feat_sd{$ft1} += ($sim - $om)*($sim - $feat_mean{$ft1});

        $feat_count{$ft2}++;
        $om = $feat_mean{$ft2};
        $feat_mean{$ft2} += (($sim - $om) / $feat_count{$ft2});
        $feat_sd{$ft2} += ($sim - $om)*($sim - $feat_mean{$ft2});

        $featpair_sim{$ft1}{$ft2} = $sim; } }

foreach (keys %feat_sd) {
    $feat_sd{$_} = sqrt($feat_sd{$_} / ($feat_count{$_}-1));
    if($feat_sd{$_} == 0) { die "\nsd is zero for $_ !\n\n"; } }


$time = runtime(); print "\n$time: normalizing and printing results ...";
# Normalizing similarity scores and printing results
open OTAB, ">$otab";
print OTAB "#feat1.id\tfeat1.desc\tfeat2.id\tfeat2.desc\tsim\tzscore\ttrues\tident\n";

my %selc_fp = ();
foreach my $ft1 (@feat1a) {
    print OMAT "$ft1\t$feat1{$ft1}";
    (my $tft1 = $ft1) =~ s/^sp[12]\.//g;

    foreach my $ft2 (@feat2a) {
        (my $tft2 = $ft2) =~ s/^sp[12]\.//g;
        
        my $s = $featpair_sim{$ft1}{$ft2};
        my $z1 = ($s - $feat_mean{$ft1}) / $feat_sd{$ft1};
        my $z2 = ($s - $feat_mean{$ft2}) / $feat_sd{$ft2};
        my $z = ($z1 + $z2) / sqrt(2);
        print OMAT "\t$z";

        print OTAB "$ft1\t$feat1{$ft1}\t$ft2\t$feat2{$ft2}";
        printf OTAB "\t%.3f\t%.3f", ($s, $z);

        my $e = join '__', sort($tft1, $tft2);
        if(exists $feat_trues{$e}) { print OTAB "\t$feat_trues{$e}"; }
        else { print OTAB "\t0"; }

        if($tft1 eq $tft2) { print OTAB "\t1\n"; }
        else { print OTAB "\t0\n"; } }
    print OMAT "\n"; }

$time = runtime(); print "\n$time: DONE\n\n";



# Index list
sub index_list {
    my $ifile = shift;
    my $itag = shift;
    my %idx = ();

    open LIST, "$ifile" or die "Can't open $ifile!";
    while (<LIST>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $idx{$itag.'.'.$p[0]} = $p[1]; }
    close LIST;

    return %idx; }


# Read-in feature-gene scores
sub get_gscore {
    my $imat = shift;
    my $itag = shift;
    my %genes = my %gscore = ();

    my $l = 0; my @col = ();
    open MAT, "$imat" or die "Can't open $imat!";
    while (<MAT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $l++; if($l==1) {
            @col = map { $_ = $itag.'.'.$_ } @p;
            next; }

        unless(exists $gene_fam->{$p[0]}) { next; }
        $genes{$p[0]}++;

        for(my $j=1; $j<=$#p; $j++) {
            $gscore{$col[$j]}{$p[0]} = $p[$j]; } }
    close MAT;

    return (\%gscore, \%genes); }


# Assigns to family to genes
sub gene2fam {
    my $ortho = shift;
    my %gfam = my %famg = ();
    my %dupg = ();

    open GFAM, "$ortho" or die "Can't open $ortho!";
    while (<GFAM>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        my $fam = shift @p;
        foreach my $spg (@p) {
            (my $g = $spg) =~ s/^[a-zA-Z]*\|//g;
            $dupg{$g}++;
            $gfam{$g} = $fam;
            $famg{$fam}{$g}++; } }
    close GFAM;

    foreach my $g (keys %gfam) {
        if($dupg{$g}>1) { delete $gfam{$g}; } }

    foreach my $f (keys %famg) {
        foreach my $g (keys %{$famg{$f}}) {
            if($dupg{$g} > 1) {
                delete $famg{$f}{$g}; } } }

    return (\%gfam, \%famg); }


# Convert gene-scores to family-scores
sub get_fscore {
    my $gscore = shift;
    my %fscore = ();
    my %fgene = ();

    my $nfam = my $ng = 0;
    foreach my $g (keys %{$gscore}) {
        my $f = $gene_fam->{$g};
        unless(exists $selc_fam{$f}) { next; }
        $ng++;
        if(exists $fscore{$f}) {
            if($gscore->{$g} > $fscore{$f}) {
                $fgene{$f} = $g;
                $fscore{$f} = $gscore->{$g}; } }
        else {
            $nfam++;
            $fgene{$f} = $g;
            $fscore{$f} = $gscore->{$g}; } }

    return (\%fscore, \%fgene); }


# Rank-transform
sub rank {
    my $href = shift;
    my @asort = sort {$a <=> $b} values %$href;
    
    my $i = 1; my %hsort = my %rep = ();
    foreach (@asort) {
        $rep{$_}++;
        if(exists $hsort{$_}) {
            $hsort{$_} = $hsort{$_} + (($i - $hsort{$_}) / $rep{$_}); }
        else { $hsort{$_} = $i; }
        $i++; }

    my %rank = ();
    foreach (keys %$href) {
        $rank{$_} = $hsort{$href->{$_}}; }
    return \%rank; }


# Calculate mean and std-deviation
sub mean_sd {
    my $href = shift;

    my $om = my $m = my $s = my $k = 0;
    foreach (keys %$href) {
        $k++;
        my $v = $href->{$_};
        $m = $om + (($v - $om) / $k);
        $s = $s + ($v - $om)*($v - $m);
        $om = $m; }

    $s = sqrt($s / ($k-1));

    return ($m, $s); }


# Calculate Spearman Rank Correlation
sub cor {
    my $rankx = shift; my $mx = shift; my $sx = shift;
    my $ranky = shift; my $my = shift; my $sy = shift;

    my %combr = ();
    my $N = scalar keys %$rankx;
    foreach my $g (keys %$rankx) {
        $combr{$g} = ($N - $rankx->{$g} + 1) * ($N - $ranky->{$g} + 1) / 2; }

    my $n = my $cor = 0;
    foreach my $g (sort {$combr{$a} <=> $combr{$b}} keys %combr) {
        $n++;
        $cor += (($rankx->{$g} - $mx) / $sx) * (($ranky->{$g} - $my) / $sy);
        if($itopn and ($n == $itopn)) { last; } }
    $cor /= ($n-1);

    return $cor; }


# Calculate Rank-biased Overlap
sub rbo {
    my $hrefx = shift;
    my $hrefy = shift;
    my $p = shift;

    my @rankx = sort {$hrefx->{$b} <=> $hrefx->{$a}} keys %$hrefx;
    my @ranky = sort {$hrefy->{$b} <=> $hrefy->{$a}} keys %$hrefy;

    my %rx = my %ry = ();
    my $rbo = 0; # my $p = 0.9;

    for(my $i=0; $i<=$#rankx; $i++) {
        $rx{$rankx[$i]}++;
        $ry{$ranky[$i]}++;

        my $c = grep { exists $rx{$_} } keys %ry;
        
        $rbo += ($c / ($i+1)) * ($p ** $i); }

    $rbo *= (1 - $p);
    return $rbo; }



__END__

=head1

Compare gene matrices from two species to return highly overlapping columns.

=head1 USAGE

compare_cross-species_matrices.pl [--imat1 INPUT_MAT1] [--imat2 INPUT_MAT2]
[--iortho ORTHOLOGY_FILE] [--imes SIMILARITY_MEASURE] [--OTAB OUTPUT_TABLE]
[--OMAT OUTPUT_MATRIX] [--help]

=head1 DESCRIPTION

This script takes in two matrices of genes (rows) x features (columns) in two
different species and compares all pairs of features at the level of
gene-families and returns their similarity scores.

=head1 ARGUMENTS

=over 12

=item C<--imat1>

Input gene x feature matrix from speicies 1.

=item C<--imat2>

Input gene x feature matrix from speicies 2.

=item C<--iortho>

Orthology information in the form of gene-families, on per line.

=item C<--imes>

(Optional) Measure of similarity. 'rbo' Rank-biased Overlap; 'cor' Spearman
Rank Correlation; Default 'rbo';

=item C<--itrues>

(Optional) True similarity scores of pairs of features.

=item C<--omat>

Output matrix containing species1_columns x species2_columns, and their
similarity scores in each cell.

=item C<--otab>

Output table containing species1_column––species2_column pairs, and their
similarity scores in each line, for pairs with normalized zscore >= 3.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2014 Feb 10

=cut

