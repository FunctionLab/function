#!/usr/bin/perl

# ================================================
# Name : smooth_scores_over_network.pl
# Purpose : Iteratively smooth node scores over network
# Created : 28-02-2014
# Last Modified : Wed 31 Dec 2014 03:11:36 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);


my $hubber = '/Genomics/ogtr04/arjunk/bin/Hubber_for_smoothing';
my $dat2dab = '/Genomics/ogtr04/arjunk/bin/Dat2Dab';

my ($help, $time, $inet, $igenes, $istep0, $itopg, $irankt, $omat);
my $initr = 20; my $iminc = 0.95; my $inetinf = 1;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'inet=s' => \$inet,
        'igenes=s' => \$igenes,
        'istep0=s' => \$istep0,
         'itopg=i' => \$itopg,
          'irankt' => \$irankt,
       'inetinf=f' => \$inetinf,
         'initr=i' => \$initr,
         'iminc=f' => \$iminc,
          'omat=s' => \$omat    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %gene_scores = ();
$gene_scores{'0'} = read_scores('0');


open MAT, ">$omat"; print MAT "#gene\titr0";
my $input_scores = 'temp.node-scores.txt';
for(my $itr=1; $itr<=$initr; $itr++) {
    print MAT "\titr$itr";
    $time = runtime(); print "\n$time: $itr ...";

    print_scores($gene_scores{$itr-1});
    run_smoother($itr);
    $gene_scores{$itr} = read_scores($itr);

    my $c = cor($gene_scores{$itr}, $gene_scores{$itr-1});
    printf " %.6f", $c;

    if($c > $iminc) { last; } }


print MAT "\n";
foreach my $g (keys %{$gene_scores{'0'}}) {
    print MAT "$g";
    foreach my $i (sort {$a <=> $b} keys %gene_scores) {
        printf MAT "\t%.6g", $gene_scores{$i}{$g}; }
    print MAT "\n"; }


`rm -f $input_scores *.step*-scores.txt`;
$time = runtime(); print "\n$time: DONE\n\n";


# Subroutines
sub read_scores {
    my $nitr = shift;

    my $scores_file;
    if($nitr==0) {
        $scores_file = $istep0; }
    else {
        ($scores_file = $inet) =~ s/\.da[tb]$//g; $scores_file =~ s/^.*\///g;
        $scores_file .= '.step'.$nitr.'-scores.txt'; }

    my %scores = (); my $sum = 0;
    open CS, "$scores_file" or die "Can't open $scores_file!";
    while (<CS>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        
        my $s = $p[$#p];
        if(($nitr > 0) and ($inetinf < 1)) {
            $s = ($s*$p[$#p-1]*$inetinf) + ((1-$inetinf)*$gene_scores{$nitr-1}{$p[0]});
            $s /= (($inetinf*$p[$#p-1]) + (1-$inetinf)); }

        $scores{$p[0]} = $s;
        $sum += $s; }
    close CS;

    foreach my $g (keys %scores) {
        $scores{$g} /= $sum; }

    return \%scores; }


sub print_scores {
    my $scores = shift;

    open IN, ">$input_scores";
    print IN "#gene\tscore\n";
    if($itopg) {
        my $nprint = 0;
        foreach my $g (sort {$scores->{$b} <=> $scores->{$a}} keys %$scores) {
            if($nprint < $itopg) {
                $nprint++; print IN "$g\t1\n"; }
            else { print IN "$g\t0\n"; } } }
    else {
        foreach my $g (keys %$scores) {
            print IN "$g\t", sprintf("%.6g\n", $scores->{$g}); } }
    close IN; }


sub run_smoother {
    my $nitr = shift;
    (my $scores_file = $inet) =~ s/\.da[tb]$//g; $scores_file =~ s/^.*\///g;
    $scores_file .= '.step'.$nitr.'-scores.txt';

    `$hubber -i $inet -w $input_scores -g -1 $igenes > $scores_file`; }


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
    my $hrefx = shift;
    my $hrefy = shift;

    my $rankx = rank($hrefx);
    my $ranky = rank($hrefy);

    my ($mx, $sx) = mean_sd($rankx);
    my ($my, $sy) = mean_sd($ranky);

    my $n = my $cor = 0;
    foreach my $g (keys %$rankx) {
        $n++;
        $cor += (($rankx->{$g} - $mx) / $sx) * (($ranky->{$g} - $my) / $sy); }
    $cor /= ($n-1);

    return $cor; }



__END__

=head1

Iteratively smooths node scores over network.

=head1 USAGE

smooth_scores_over_network.pl.pl [--inet INPUT_NETWORK] [--igenes NET_GENES]
[--istep0 INITIAL_NODE-SCORES] [--omat OUTPUT_MAT] [--help]

=head1 DESCRIPTION

This script takes in a network and node-scores, and computes new node-scores
that have benn smoothed over the network itertatively.

=head1 ARGUMENTS

=over 12

=item C<--inet>

Input network in DAT/DAB format.

=item C<--igenes>

List of genes in the network.

=item C<--istep0>

Initial node-scores in 2 tab-seprated columns: <gene> <score>.

=item C<--initr>

(Optional) No. of iterations. Default 20.

=item C<--iminc>

(Optional) Min. correlation to stop smoothing. Default 0.95.

=item C<--omat>

Output matrix containing genes along the rows and iteratively updated scores
along the columns.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2014 Feb 28

=cut

