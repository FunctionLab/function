#!/usr/bin/perl

# ================================================
# Name : calc_shared-nearest-neighbors.pl
# Purpose : Calculate edge weights based on shared top k neighbors of genes
# Created : 06-01-2015
# Last Modified : Thu 26 Feb 2015 01:01:19 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use lib '/Genomics/Users/arjunk/software/lib/perl5/';
use PDF;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);


my ($help, $idat, $icut, $odat);
my $itopk = 50; my $ipval = 1E-8;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions(     'help' => \$help,
              'idat=s' => \$idat,
             'itopk=i' => \$itopk,
              'icut=f' => \$icut,
             'ipval=f' => \$ipval,
              'odat=s' => \$odat    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my $time = runtime(); print "\n$time: Parsing DAT and Indexing topk genes ...";

#my %tempn = ();
my %topk_genes = my %topk_values = my %topk_count = ();
my $one = 0;
open IDAT, "$idat" or die "Can't open $idat!";
while (<IDAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my ($g1, $g2, $s) = split '\t', $_;
    print "\n=> $g1 $g2 $s";

    if($icut and ($s < $icut)) { next; }
    $one++;
    unless($one % 100000) {
        $time = runtime(); print "\n\t$time: $one edges ..."; }

    unless(exists $topk_count{$g1}) { $topk_count{$g1} = 0; }
    if($topk_count{$g1} < $itopk) {
        $topk_count{$g1}++;
        push(@{$topk_genes{$g1}}, $g2);
        push(@{$topk_values{$g1}}, $s);
        print "\n$g1 <- push $g2 $s : $topk_count{$g1}"; }
    else {
        print "\n$g1 : topk : ";
        print join " ", @{$topk_genes{$g1}}, ":", @{$topk_values{$g1}}; }

    if($topk_count{$g1} == $itopk) {
        G1: for(my $i=($itopk-1); $i>=0; $i--) {
            print "\n$g1 : $i :";
            if($s <= ${$topk_values{$g1}}[$i]) {
                if($i < ($itopk-1)) {
                    print " splice $g2 next to $i and pop";
                    splice(@{$topk_genes{$g1}}, ($i+1), 0, $g2); pop @{$topk_genes{$g1}};
                    splice(@{$topk_values{$g1}}, ($i+1), 0, $s); pop @{$topk_values{$g1}}; }
                else { print " $g2 $s smaller than last"; }
                last G1; } } }

    unless(exists $topk_count{$g2}) { $topk_count{$g2} = 0; }
    if($topk_count{$g2} < $itopk) {
        $topk_count{$g2}++;
        push(@{$topk_genes{$g2}}, $g1);
        push(@{$topk_values{$g2}}, $s);
        print "\n$g2 <- push $g1 $s : $topk_count{$g2}"; }
    else {
        print "\n$g2 : topk : ";
        print join " ", @{$topk_genes{$g2}}, ":", @{$topk_values{$g2}}; }

    if($topk_count{$g2} == $itopk) {
        G2: for(my $i=($itopk-1); $i>=0; $i--) {
            print "\n$g2 : $i :";
            if($s <= ${$topk_values{$g2}}[$i]) {
                if($i < ($itopk-1)) {
                    print " splice $g1 next to $i and pop";
                    splice(@{$topk_genes{$g2}}, ($i+1), 0, $g1); pop @{$topk_genes{$g2}};
                    splice(@{$topk_values{$g2}}, ($i+1), 0, $s); pop @{$topk_values{$g2}}; }
                else { print " $g1 $s smaller than last"; }
                last G2; } } } }
close IDAT;

my @agenes = sort keys %topk_genes;
my $ong = scalar @agenes;
print "\n\t$ong genes; $one edges";

my %topk = ();
foreach my $g1 (@agenes) {
    print "\n$g1 ->";
    foreach my $g2 (@{$topk_genes{$g1}}) {
        print " $g2";
        $topk{$g1}{$g2}++; } }
print "\n"; exit;


$time = runtime(); print "\n$time: Calc & printing new edge weights ...";

open ODAT, ">$odat";
my $ne = 0;

for(my $i=0; $i<$#agenes; $i++) {
    my $g1 = $agenes[$i];
    my $size1 = scalar keys %{$topk{$g1}};

    for(my $j=($i+1); $j<=$#agenes; $j++) {
        my $g2 = $agenes[$j];
        my $size2 = scalar keys %{$topk{$g2}};

        $ne++;
        unless($ne % 100000) {
            $time = runtime(); print "\n\t$time: $ne edges ..."; }

        my @comg = grep { exists $topk{$g1}{$_} } keys %{$topk{$g2}};
        my $ncom = scalar @comg;

        my $pval = hypergeometric_tail($ong, $size1, $size2, $ncom);
        if($pval > $ipval) { next; }

        print ODAT "$g1\t$g2\t", sprintf("%.6g\n", (-1*log($pval)/log(10)));
    }
}

close ODAT;

$time = runtime(); print "\n$time: DONE\n\n";



__END__

=head1

Calculate new edge weights based on shared top k neighbors of genes

=head1 USAGE

calc_shared-nearest-neighbors.pl [--idat NETWORK_DAT] [--itopk TOP-K_NEIGHBORS] [--icut
EDGE_CUTOFF] [--ipval PVALUE_CUTOFF] [--odat NETWORK_DAT] [--help]

=head1 DESCRIPTION

This script takes in a network in DAT format, and calculates new edges between
nodes based on their shared k-nearest neighbors. The edge weights are equal to
the negative log of the hypergeometric p-value of overlap between the shared
neighbors of a pair of nodes. Only edges with p-value < --ipval are printed. If
provided, --icut is used consider neighbors only >= this cutoff.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input network in DAT format. This algorithm works much better if the edges in
the DAT are sorted in decreasing edge weight.

=item C<--itopk>

(Optional) Number of top edges to retain per node. Default is 50.

=item C<--icut>

(Optional) Edge weight cutoff.

=item C<--ipval>

(Optional) New edge p-value cutoff. Default is 1E-8.

=item C<--odat>

Output network in DAT format.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2015 Jan 05

=cut

