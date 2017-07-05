#!/usr/bin/perl

# ================================================
# Name : pcl2ssv.pl
# Purpose : Convert PCL profile into SSV format
# Created : 15-05-2015
# Last Modified : Sun 12 Jun 2016 02:44:15 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);

sub quant_norm;
sub gene_rank;


my ($help, $ipcl, $isam, $inzo);
my $igen = '/Genomics/ogtr04/arjunk/projects/human-sample-classification/age-group/standards/mas5_gpl570_genes.txt';
my $iref = '/Genomics/ogtr04/arjunk/data/gene-expression/gse/sample-annotaions/age-group-evaluation/data/upc-mas5-comparison/mas5_genome-ref.txt';

pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(     'help' => \$help,
              'ipcl=s' => \$ipcl,
              'isam=s' => \$isam,
              'igen=s' => \$igen,
              'iref=s' => \$iref,
                'inzo' => \$inzo    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %genes = my @agenes = ();
open GEN, "$igen" or die "Can't open $igen!";
while (<GEN>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $genes{$_}++; push(@agenes, $_); }
close GEN;


my @refprof = ();
open REF, "$iref" or die "Can't open $iref!";
chomp(@refprof = <REF>); close REF;
my $ref_len = scalar @refprof;

@refprof = sort {$a <=> $b} @refprof;

unless($inzo) {
    my $lwr = $refprof[int(($ref_len*0.0005) + 0.5) - 1];
    my $hgr = $refprof[int(($ref_len*0.9995) + 0.5) - 1];

    foreach (@refprof) {
        if($_ < $lwr) { $_ = $lwr; }
        elsif($_ > $hgr) { $_ = $hgr; }
        $_ = ($_ - $lwr)/($hgr - $lwr); } }

my $ref_median = median(\@refprof);


my %sex_sample = ();
if($isam) {
    open SAM, "$isam" or die "Can't open $isam!";
    while (<SAM>) {
        if($_ =~ /^#/) { next; }
        if($_ =~ /male\.score/) { next; }
        chomp($_); my @p = split '\t', $_;

        my $xpre;
        if($p[$#p-1] > $p[$#p]) { $xpre = 'female'; }
        else { $xpre = 'male'; }

        $sex_sample{$xpre}{$p[0]}++; }
    close SAM; }


my @asamples = my %sample_gene_exp = (); my $l = 0;
open PCL, "$ipcl" or die "Can't open $ipcl!";
while (<PCL>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    $l++;
    if($l == 1) {
        @asamples = @p; next; }

    unless(exists $genes{$p[0]}) { next; }
    for(my $j=1; $j<=$#p; $j++) {
        unless(looks_like_number($p[$j])) { next; }
        $sample_gene_exp{$asamples[$j]}{$p[0]} = $p[$j]; } }
close PCL;


(my $ossv = $ipcl) =~ s/^.*\///g; $ossv =~ s/\.pcl$/\.ssv/g;
(my $omssv = $ipcl) =~ s/^.*\///g; $omssv =~ s/\.pcl$/\.male.ssv/g;
(my $ofssv = $ipcl) =~ s/^.*\///g; $ofssv =~ s/\.pcl$/\.female.ssv/g;
#print "\n\n$ossv\n$omssv\n$ofssv\n"; exit;

if($isam) {
    foreach my $x (keys %sex_sample) {
        if($x =~ /^m/) { open MSSV, ">$omssv"; }
        elsif($x =~ /^f/) { open FSSV, ">$ofssv"; } } }
else {
    open SSV, ">$ossv"; }

shift @asamples; my $nsample = 0;
foreach my $s (@asamples) {
    $nsample++;
    #my $ossv = "s$nsample.ssv";
    #open SSV, ">$ossv";

    print "\n$nsample\t$s ...";
    my $feat = quant_norm($sample_gene_exp{$s}, \@refprof, $ref_median, \@agenes);
    #close SSV;

    if($isam) {
        if(exists $sex_sample{'male'}{$s}) {
            print MSSV "0 ", (join ' ', @{$feat}), "\n"; }
        elsif(exists $sex_sample{'female'}{$s}) {
            print FSSV "0 ", (join ' ', @{$feat}), "\n"; } }
    else {
        print SSV "0 ", (join ' ', @{$feat}), "\n"; } }

if($isam) {
    if(-e $omssv) { close MSSV; }
    if(-e $ofssv) { close FSSV; } }
else {
    close SSV; }
print "\n\n";


# Get gene ranks
sub gene_rank {
    my $gref = shift;

    my @asort = sort {$a <=> $b} values %$gref;
    my $N = scalar @asort;
    
    my $i = 1; my %hsort = my %rep = ();
    foreach (@asort) {
        $rep{$_}++;
        if(exists $hsort{$_}) {
            $hsort{$_} = $hsort{$_} + (($i - $hsort{$_}) / $rep{$_}); }
        else { $hsort{$_} = $i; }
        $i++; }

    my %rank = (); my $minr = $N; my $maxr = 1;
    foreach my $g (keys %$gref) {
        $rank{$g} = $hsort{$gref->{$g}};
        if($rank{$g} < $minr) { $minr = $rank{$g}; }
        if($rank{$g} > $maxr) { $maxr = $rank{$g}; } }

    return (\%rank, $minr, $maxr); }


# Quantile-transform
sub quant_norm {
    my $gref = shift; # input gene exp
    my $rref = shift; # ref. profile
    my $rmed = shift; # median of ref. profile
    my $aref = shift; # array of genes

    my $Nref = scalar @$rref;
    my ($grnk, $minr, $maxr) = gene_rank($gref);

    my @qnorm = (); my $j = 0;
    foreach (@$aref) {
        $j++;
        if(exists $grnk->{$_}) {
            my $r = $grnk->{$_};
            my $ref_idx = int(($Nref - 1)*($r - $minr)/($maxr - $minr) + 0.5);
            push(@qnorm, ($j.':'.sprintf("%.6f", ${$rref}[$ref_idx]))); }
        else {
            push(@qnorm, ($j.':'.sprintf("%.6f", $rmed))); } }

    return \@qnorm; }


# Median
sub median {
    my $aref = shift;
    my @ary = sort {$a <=> $b} grep {!/^NA$/} @{$aref};
    if(@ary % 2) { return $ary[@ary/2]; }
    else { return (($ary[(@ary/2)-1] + $ary[@ary/2]) / 2); } }




__END__

=head1

Convert PCL to SSV for prediction.

=head1 USAGE

pcl2ssv.pl [--ipcl INPUT_PCL] [--igen GENE_FEATURES] [--iref REFERENCE_PROFILE] [--inzo] [--help]

=head1 DESCRIPTION

This script takes in a gene expression matrix in tabular format and converts it into
a SSV format containing samples as examples/instances along the rows and
gene-expression values as features, to be used for prediction using liblinear.
The output is stored in the <input>.ssv file.

=head1 ARGUMENTS

=over 12

=item C<--ipcl>

Input gene-expression matrix with the first row containing sample identifiers,
first column containing Entrez gene identifiers, and the rest of the matrix
containing gene-expression values. NA's can be used to denote missing values.

=item C<--isam>

(Optional) Sample information in tabular form that is an output of the
predict-sex script or similar with sex of each sample in the penultimate column.

=item C<--igen>

(Optional) Complete list of genes to be used as the feature vector.
Default: /Genomics/ogtr04/arjunk/projects/human-sample-classification/age-group/standards/mas5_gpl570_genes.txt

=item C<--iref>

(Optional) Reference profile to be used as target for quantile transformation.
Default: /Genomics/ogtr04/arjunk/data/gene-expression/gse/sample-annotaions/age-group-evaluation/data/upc-mas5-comparison/mas5_genome-ref.txt

=item C<--inzo>

(Optional) Indicate that the quantile-transformed values need not be scaled
between 0 and 1. If not provided, values are linearly 0-1 scaled.

=item C<--help>

Prints this help message.

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

