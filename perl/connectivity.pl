#!/usr/bin/perl

# ================================================
# Name : connectivity.pl
# Purpose : Gene function prediction using local network connectivity
# Created : 18-03-2013
# Last Modified : Mon 03 Mar 2014 11:31:24 AM EST
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
sub evaluate;

my $hubber = '/r04/arjunk/bin/Hubber';

my ($help, $inet, $iglab, $igfilt, $ipred, $otab, $time, @p);
my $ineg = 1; my $icvk = 5; my $icvn = 10;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions( 'help' => \$help,
          'inet=s' => \$inet,
         'iglab=s' => \$iglab,
        'igfilt=s' => \$igfilt,
          'ineg=s' => \$ineg,
          'icvk=i' => \$icvk,
          'icvn=i' => \$icvn,
           'ipred' => \$ipred,
          'otab=s' => \$otab ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);
pod2usage( -exitstatus => 2, -verbose => 2 ) unless ($inet and $iglab and $otab);

(my $t1 = $iglab) =~ s/^.*\///g; $t1 =~ s/\.(genes|txt)//g;
(my $t2 = $inet) =~ s/^.*\///g; $t2 =~ s/\.dab$//g;


# Parse selected genes
my %sel_allg = ();
if($igfilt) {
    %sel_allg = get_selec($igfilt); }


# Parse gene labels
my %gene_truel = my %truel_gene = ();
open GL, "$iglab" or die "Can't open $iglab!";
while (<GL>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    if($igfilt) { unless(exists $sel_allg{$p[0]}) { next; } }

    $gene_truel{$p[0]} = $p[1];
    if($p[1]==0) { next; }
    $truel_gene{$p[1]}{$p[0]}++; }
close GL;


if($ipred) {
    $time = runtime(); print "\n$time: Printing standard ...\n";
    my $dir = $t1.'.'.$t2.'.neg'.$ineg; `mkdir -p $dir`;

    foreach my $lab (keys %truel_gene) {
        if(($ineg==0) and ($lab==-1)) { next; }
        my $glist = $dir.'/'.$lab;
        open GH, ">$glist";
        foreach my $g (keys %{$truel_gene{$lab}}) {
            print GH "$g\n"; }
        close GH; }


    $time = runtime(); print "\n$time: Running hubber ...\n";
    my $out = $t1.'.'.$t2.'.neg'.$ineg.'.hubber.out';
    `$hubber -i $inet -g -1 $dir/* > $out`;


    $time = runtime(); print "\n$time: Recording predictions ...\n";
    my %gene_lab_preds = (); my @genes;
    open PS, "$out" or die "Can't open $out!";
    while (<PS>) {
        chomp($_); @p = split '\t', $_;
        if($_ =~ /^Function/) {
            @genes = @p; next; }

        for(my $j=1; $j<=$#p; $j++) {
            my $g = $genes[$j];
            $gene_lab_preds{$g}{$p[0]} = abs($p[$j]); } }
    close PS;

    my %gene_preds = ();
    foreach my $g (keys %gene_lab_preds) {
        my $preds;
        if($ineg) {
            $preds = $gene_lab_preds{$g}{'1'} -
                $gene_lab_preds{$g}{'-1'}; }
        else {
            $preds = $gene_lab_preds{$g}{'1'}; }

        $gene_preds{$g} = $preds; }


    $time = runtime(); print "\n$time: Printing predictions ...\n";
    open my $TAB, ">$otab";
    print $TAB "#$t1.$t2.neg$ineg\n";
    foreach my $g (sort {$gene_preds{$b} <=> $gene_preds{$a}} keys %gene_preds) {
        print $TAB "$g\t$gene_truel{$g}\t$gene_preds{$g}\n"; }

    $time = runtime(); print "\n$time: DONE\n\n";
    exit; }


$time = runtime(); print "\n$time: Splitting standard for CV ...\n";
# Splitting standard for CV
my $dir = $t1.'.'.$t2.'.neg'.$ineg.'.'.$icvn.'-'.$icvk.'-cv'; `mkdir -p $dir`;
print "\tsaving sets in $dir/\n";

my %gene_cvi_train = my %gene_cvi_eval = ();
foreach my $lab (keys %truel_gene) {
    my %cvi_eval = get_cvint($truel_gene{$lab}, $icvn, $icvk);
    foreach my $cvi (keys %cvi_eval) {
        my $glist = $dir.'/'.$cvi.'_'.$lab;
        if(($ineg) or (($ineg==0) and ($lab==1))) {
            open GH, ">$glist";

            foreach my $g (keys %{$truel_gene{$lab}}) {
                if(exists $cvi_eval{$cvi}{$g}) {
                    $gene_cvi_eval{$g}{$cvi}++;
                    next; }
                $gene_cvi_train{$g}{$cvi}++;
                print GH "$g\n"; }

            close GH; } } }


$time = runtime(); print "\n$time: Running hubber ...\n";
# Running Hubber
my $out = $t1.'.'.$t2.'.neg'.$ineg.'.'.$icvn.'-'.$icvk.'-cv.hubber.out';
print "\tsaving results to $out\n";
`$hubber -i $inet -g -1 $dir/* > $out`;


$time = runtime(); print "\n$time: Recording predictions ...\n";
# Recording predicted gene scores for each function
my %gene_cvi_lab_preds = (); my @genes;
open CV, "$out" or die "Can't open $out!";
while (<CV>) {
    chomp($_); @p = split '\t', $_;
    if($_ =~ /^Function/) {
        @genes = @p; next; }

    my ($cvi, $lab) = split '_', $p[0];
    
    for(my $j=1; $j<=$#p; $j++) {
        my $g = $genes[$j];
        if(exists $gene_cvi_train{$g}{$cvi}) { next; }
        $gene_cvi_lab_preds{$g}{$cvi}{$lab} = $p[$j]; } }
close CV;

my %gene_preds = my %gene_predc = ();
foreach my $g (keys %gene_cvi_lab_preds) {
    $gene_preds{$g} = $gene_predc{$g} = 0;

    foreach my $cvi (keys %{$gene_cvi_lab_preds{$g}}) {
        my $preds;
        if($ineg) {
            $preds = $gene_cvi_lab_preds{$g}{$cvi}{'1'} -
                $gene_cvi_lab_preds{$g}{$cvi}{'-1'}; }
        else {
            $preds = $gene_cvi_lab_preds{$g}{$cvi}{'1'}; }

        $gene_predc{$g}++;
        $gene_preds{$g} += ($preds - $gene_preds{$g})/$gene_predc{$g}; } }


$time = runtime(); print "\n$time: Evaluating and Printing results ...\n";
open my $TAB, ">$otab";
print $TAB "#$t1.$t2.neg$ineg\n";
foreach my $g (sort {$gene_preds{$b} <=> $gene_preds{$a}} keys %gene_preds) {
    print $TAB "$g\t$gene_truel{$g}\t$gene_preds{$g}\t$gene_predc{$g}\n"; }


evaluate(\%gene_truel, \%gene_preds, $TAB);
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

# Given indices and size for a class, assigns indices to cv-fold
sub get_cvint {
    my $egs = shift;
    my $cvn = shift;
    my $cvk = shift;

    my %cvint = ();

    for(my $n=0; $n<$cvn; $n++) {
        my @rnd = shuffle( keys %$egs );
        my $i = 0;

        foreach (@rnd) {
            my $cvi = $n.'-'.($i % $cvk);
            $cvint{$cvi}{$_}++;
            $i++; } }

    return %cvint; }

# Randomize array
sub fy_shuffle {
    my $array = shift;
    my $i;

    srand(10*$icvn);
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i]; } }

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

    print $FH "#P\t$P\n#N\t$N\n";
    print $FH "#PRIOR\t", sprintf("%.3f\n", ($P/($P+$N)));
    print $FH "#AUC\t", sprintf("%.3f\n", $auc);
    print $FH "#P10R\t", sprintf("%.3f\n", $p10r);
    print $FH "#P20R\t", sprintf("%.3f\n", $p20r);
    print $FH "#P50R\t", sprintf("%.3f\n", $p50r);
    print $FH "#AUPRC\t", sprintf("%.3f\n", $auprc);
    print $FH "#WAUPRC\t", sprintf("%.3f\n", $wauprc);
}


__END__

=head1

Run network-based gene function prediction using local network connectivity.

=head1 USAGE

connectivity.pl [--inet NETWORK_DAB] [--iglab GENE_LABELS] [--igfilt
GENE_FILTER] [--ineg NEG_FLAG] [--icvk CROSS-VALIDATION_FOLDS] [--icvn
CROSS-VALIDATION_REPEATS] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in a network and gene labels and runs network-based
gene-function-prediction using local connectivity. This method can use negatives
along with positives during prediction in addition to just evaluation.

=head1 ARGUMENTS

=over 12

=item C<--inet>

Network file in DAB format.

=item C<--iglab>

List of genes and their labels in a tab-delimited format: <gene> <-1/0/1>
Genes labelled '0' will be used just for prediction.

=item C<--igfilt>

(Optional) List of genes to include.

=item C<--ineg>

(Optional) Flag to determine if negatives will be used in prediction. '0' if
–ves are to be used only for evaluation. '1' if they are to be used for
prediction too. Default '1'.

=item C<--icvk>

(Optional) Cross-validation fold. Default 5.

=item C<--icvn>

(Optional) Cross-validation repeats. Default 10.

=item C<--ipred>

(Optional) Flag to only do prediction based on labelled examples without any CV.

=item C<--otab>

Output table containing prediction and evaluation statistics.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 Mar 03

=cut

