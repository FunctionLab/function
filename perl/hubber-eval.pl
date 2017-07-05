#!/usr/bin/perl

# ================================================
# Name : hubber-predict.pl
# Purpose : Evaluate gene function prediction using Hubber
# Created : 18-03-2013
# Last Modified : Wed 31 Dec 2014 12:44:33 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);
use List::Util qw(shuffle);

sub get_selec;
sub parse_gmt;
sub get_cvint;
sub dcheck;

my ($help, $inet, $igmt, $igsf, $islim, $igenef, $otab, $time, @p);
my $imings = 5; my $imaxgs = 100; my $icvk = 5; my $icvn = 10;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions( 'help' => \$help,
          'inet=s' => \$inet,
          'igmt=s' => \$igmt,
          'igsf=s' => \$igsf,
           'islim' => \$islim,
        'igenef=s' => \$igenef,
        'imings=s' => \$imings,
        'imaxgs=s' => \$imaxgs,
          'icvk=i' => \$icvk,
          'icvn=i' => \$icvn,
          'otab=s' => \$otab ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);
pod2usage( -exitstatus => 2, -verbose => 2 ) unless ($inet and $igmt and $otab);


$time = runtime(); print "\n$time: Parsing genes and genesets ...";

# Parse selected genes
my %sel_genes = ();
if($igenef) {
    %sel_genes = get_selec($igenef); }

# Parse GMT
my @ref = parse_gmt($igmt, \%sel_genes);
my %gs_posg = %{$ref[0]};
my %gs_desc = %{$ref[1]};
my %gs_size = %{$ref[2]};
my %gs_posg_cvint = %{$ref[3]};

$time = runtime(); print "\n$time: Obtaining +ves and –ves ...";

# Assign true gene-labels and record all genes
my %gs_gene_truel = my %all_genes = ();
foreach my $gs (keys %gs_posg) {
    foreach my $g (keys %{$gs_posg{$gs}}) {
        $gs_gene_truel{$gs}{$g} = 1;
        $all_genes{$g}++; } }

# Record slim-term to gene mapping
my %gs_slim = my %slim_genes = ();
if($islim) {
    %gs_slim = %{$ref[4]};
    foreach my $gs (keys %gs_slim) {
        foreach my $s (keys %{$gs_slim{$gs}}) {
            foreach my $g (keys %{$gs_posg{$gs}}) {
                $slim_genes{$s}{$g}++; } } } }

# Parse selected genesets
my %sel_gs = ();
if($igsf) {
    %sel_gs = get_selec($igsf); }

# Get positives and negatives for each geneset
my %gs_negg = my %gs_negg_cvint = ();
my %genes_to_avoid;
foreach my $gs (keys %gs_posg) {
    if($igsf) { unless(exists $sel_gs{$gs}) { next; } }
    %genes_to_avoid = ();

    if($islim) {
        foreach my $s (keys %{$gs_slim{$gs}}) {
            foreach my $g (keys %{$slim_genes{$s}}) {
                $genes_to_avoid{$g}++; } } }

    foreach my $g (keys %all_genes) {
        if(exists $genes_to_avoid{$g}) { next; }
        if(exists $gs_posg{$gs}{$g}) { next; }
        $gs_gene_truel{$gs}{$g} = -1;
        $gs_negg{$gs}{$g}++; }

    $gs_negg_cvint{$gs} = get_cvint($gs_negg{$gs}, $icvk); }

$time = runtime(); print "\n$time: Cross-validation ...";

# CV
(my $dir = $inet) =~ s/^.*\///g; $dir =~ s/\.dab$//g;
(my $mat = $inet) =~ s/^.*\///g; $mat =~ s/\.dab$//g;
my ($cvdir, $cvmat);

my (@genes, $score); my %gs_gene_score = ();

for(my $k=0; $k<=$icvk; $k++) {
    $cvdir = $dir.'_cv'.$k;
    `mkdir -p $cvdir`;
    $time = runtime(); print "\n\t$time: $cvdir ...";

    # Printing contexts
    foreach my $gs (keys %sel_gs) {
        open GH, ">$cvdir/$gs";
        foreach my $g (keys %{$gs_posg{$gs}}) {
            if(exists $gs_posg_cvint{$gs}{$k}{$g}) { next; }
            print GH "$g\n"; }
        close GH; }

    $time = runtime(); print "\n\t\t$time: Running hubber ...";

    # Running Hubber
    $cvmat = $mat.'_cv'.$k.'.hubber.out';
    `Hubber -i $inet -g -1 $cvdir/* > $cvmat`;

    $time = runtime(); print "\t\t$time: Recording predictions ...";

    # Recording in predicted gene scores for each function
    open CV, "$cvmat" or die "Can't open $cvmat!";
    while (<CV>) {
        chomp($_); @p = split '\t', $_;
        if($_ =~ /^Function/) {
            @genes = @p; next; }
        
        my $gs = $p[0];
        for(my $j=1; $j<=$#p; $j++) {
            my $g = $genes[$j];
            if((exists $gs_posg_cvint{$gs}{$k}{$g})
                    or (exists $gs_negg_cvint{$gs}{$k}{$g})) {
                if($p[$j] < 0) { die "unexp -ve! $gs $g\n"; }
                if(exists $gs_gene_score{$gs}{$g}) { die "present! $gs $g\n"; }
                $gs_gene_score{$gs}{$g} = $p[$j]; } } }
    close CV;

    `rm -rf ./$cvdir`;
    `rm -f $cvmat`; }

$time = runtime(); print "\n$time: Evaluating and Printing results ...";

# Print evaluation results
open my $TAB, ">$otab";
print $TAB "Func\tDesc\t#P\t#N\tAUC\tP10R\tP20R\tP50R\tAUPRC\tWAUPRC\n";
foreach my $gs (keys %gs_gene_score) {
    dcheck($gs, $gs_gene_truel{$gs}, $gs_gene_score{$gs}, $TAB); }
close $TAB;

$time = runtime(); print "\n$time: DONE\n\n";


# Subroutines
sub get_selec {
    my $f = shift;
    my %sel = (); my @q;

    open SF, "$f" or die "Can't open $f!";
    while (<SF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @q = split '\t', $_;
        $sel{$q[0]}++; }
    close SF;

    return %sel; }

sub parse_gmt {
    my $gmt = shift;
    my $bgr = shift;
    my ($gs, $s, @q, @r);

    my %genes = my %desc = my %size = my %cvint = ();
    my %slim = ();

    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @q = split '\t', $_;
        
        $gs = shift @q;
        ($desc{$gs} = shift @q) =~ s/ \([0-9]*\)$//g;

        if($islim) {
            @r = split ' ', $desc{$gs};
            $s = pop @r;
            $desc{$gs} = join ' ', @r;

            if($s =~ /\|/) {
                @r = split '\|', $s;
                foreach my $a (@r) {
                    $slim{$gs}{$a}++; } }
            else { $slim{$gs}{$s}++; } }
        
        $size{$gs} = 0;
        foreach my $g (@q) {
            if($igenef) {
                unless(exists $bgr->{$g}) { next; } }
            $genes{$gs}{$g}++;
            $size{$gs}++; }

        if($size{$gs} == 0) {
            delete $desc{$gs};
            delete $size{$gs};
            if($islim) { delete $slim{$gs}; }
            next; }

        if(($size{$gs} < $imings) or ($size{$gs} > $imaxgs)) { next; }
        $cvint{$gs} = get_cvint($genes{$gs}, $icvk);
    }
    close GMT;

    return (\%genes, \%desc, \%size, \%cvint, \%slim); }

# Given indices and size for a class, assigns indices to cv-fold
sub get_cvint {
    my $dat = shift;
    my $cvk = shift;

    # Randomize the indices
    my @rnd = shuffle( keys %$dat );

    my %cvi = (); my $i = 0;

    foreach (@rnd) {
        $cvi{$i++ % $cvk}{$_}++; }

    return \%cvi; }

# Randomize array
sub fy_shuffle {
    my $array = shift;
    my $i;

    srand(10*$icvn);
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i]; } }

sub dcheck {
    my $funci = shift;
    my $truel = shift;
    my $preds = shift;
    my $FH = shift;

    my $P = my $N = 0;
    foreach (keys %$truel) {
        if($truel->{$_} == 1) { $P++; }
        else { $N++; } }

    my $auc = my $auprc = my $wauprc = my $p10r = my $p20r = my $p50r = 0;
    my $tp = my $fp = my $tn = my $fn = 0;
    my $in10 = my $in20 = my $in50 = 1;
    my $pr = my $rc = my $w = my $sumw = 0;
    
    foreach (sort {$preds->{$b} <=> $preds->{$a}} keys %{$preds}) {
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

            $w = (-1)*log($rc) / log(2);
            $sumw += $w;
            $wauprc += $w*$pr;
        }
        else {
            $fp++;
            $pr = $tp / ($tp + $fp);
            $auc += $tp; }

        $tn = $N - $fp; $fn = $P - $tp; }

    if(($tp == 0) or ($fp == 0)) {
        print "warning: Too few +ve true labels or -ve true labels\n";
        $auc = $auprc = $wauprc = 0; }
    else {
        $auc = $auc / $tp / $fp;
        $auprc = $auprc / $tp;
        $wauprc = $wauprc / $sumw; }

    print $FH "$funci\t$gs_desc{$funci}\t$P\t$N\t", sprintf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
        $auc, $p10r, $p20r, $p50r, $auprc, $wauprc), "\n"; }


__END__

=head1

Run evaluation of network-based gene function prediction using Hubber.

=head1 USAGE

hubber-predict.pl [--inet NETWORK_DAB] [--igmt GENESETS_GMT] [--igsf
GENESET_FILTER] [--islim] [--igenef GENE_FILTER] [--imings MIN_GENESET_SIZE]
[--imaxgs MAX_GENESET_SIZE] [--icvk CROSS-VALIDATION_FOLDS] [--icvn
CROSS-VALIDATION_REPEATS] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in a network and a collection of genesets and runs
evaluation of gene-function-prediction using Hubber.

=head1 ARGUMENTS

=over 12

=item C<--inet>

Network file in DAB format.

=item C<--igmt>

Geneset collection in GMT format.

=item C<--igsf>

(Optional) List of genesets in evaluate (one per line).

=item C<--islim>

(Optional) Set this option when the GMT file contains slim terms appended to
geneset's description: <gs_id><\t><gs_desc> <slim1|slim2>
(<no.genes>)<\t><gene1><\t><gene2> ...

=item C<--igenef>

(Optional) List of genes to include.

=item C<--imings>

(Optional) Genesets with less than --imings genes will be ignored. Default 5.

=item C<--imaxgs>

(Optional) Genesets with more than --imaxgs genes will be ignored. Default 100.

=item C<--icvk>

(Optional) Cross-validation fold. Default 5.

=item C<--icvn>

(Optional) Cross-validation repeats. Default 10.

=item C<--otab>

Output table containing evaluation statistics.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 Mar 03

=cut

