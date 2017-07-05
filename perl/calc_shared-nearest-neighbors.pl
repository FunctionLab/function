#!/usr/bin/perl

# ================================================
# Name : calc_shared-nearest-neighbors.pl
# Purpose : Calculate edge weights based on shared top k neighbors of genes
# Created : 06-01-2015
# Last Modified : Mon 02 Mar 2015 03:44:52 PM EST
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
my $itopk = 50; my $ipval = 1E-6;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions(     'help' => \$help,
              'idat=s' => \$idat,
             'itopk=i' => \$itopk,
              'icut=f' => \$icut,
             'ipval=f' => \$ipval,
              'odat=s' => \$odat    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my $time = runtime(); print "\n$time: Parsing DAT and Indexing topk genes ...";

my %topk_genes = my %topk_count = ();
my $one = 0;
open IDAT, "$idat" or die "Can't open $idat!";
while (<IDAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my ($g1, $g2, $s) = split '\t', $_;

    if($icut and ($s < $icut)) { next; }
    $one++;

    unless(exists $topk_count{$g1}) { $topk_count{$g1} = 0; }
    if($topk_count{$g1} < $itopk) {
        $topk_count{$g1}++;
        push(@{$topk_genes{$g1}}, $g2); }

    unless(exists $topk_count{$g2}) { $topk_count{$g2} = 0; }
    if($topk_count{$g2} < $itopk) {
        $topk_count{$g2}++;
        push(@{$topk_genes{$g2}}, $g1); }
}
close IDAT;

my @tempg = sort {$topk_count{$b} <=> $topk_count{$a}} keys %topk_genes;

my %topk = my @agenes = my %all_topk = (); my $ong_ltk = 0;
foreach my $g1 (@tempg) {
    if($topk_count{$g1} < 4) { next; }
    if($topk_count{$g1} < $itopk) { $ong_lek++; }
    push(@agenes, $g1);
    foreach my $g2 (@{$topk_genes{$g1}}) {
        $all_topk{$g2}++;
        $topk{$g1}{$g2}++; } }
my $ong = scalar @agenes;
print "\n\t$ong genes; $one edges";
print "\n\tnp. genes with topk < $itopk: $ong_ltk";


$time = runtime(); print "\n$time: Calc & printing new edge weights ...";

open ODAT, ">$odat";

if($odat =~ /\.dat$/) {
    #my $ne = 0;
    for(my $i=0; $i<$#agenes; $i++) {
        my $g1 = $agenes[$i];
        my $size1 = scalar keys %{$topk{$g1}};

        for(my $j=($i+1); $j<=$#agenes; $j++) {
            my $g2 = $agenes[$j];
            my $size2 = scalar keys %{$topk{$g2}};

            #$ne++;
            #unless($ne % 1000000) {
            #    $time = runtime(); print "\n\t$time: $ne edges ..."; }

            my @comg = grep { exists $topk{$g1}{$_} } keys %{$topk{$g2}};
            my $ncom = scalar @comg;
            if($ncom < 3) { next; }

            my $pval = hypergeometric_tail($ong, $size1, $size2, $ncom);
            if($pval > $ipval) { next; }

            print ODAT "$g1\t$g2\t", sprintf("%.6f\n", (-1*log($pval)/log(10))); } } }

if($odat =~ /\.pcl$/) {
    my @tgenes = sort keys %all_topk;
    foreach (@tgenes) {
        print ODAT "\t$_"; }
    print ODAT "\n";

    foreach my $g (@agenes) {
        print ODAT "$g";

        foreach my $h (@tgenes) {
            my $v = 0;
            if(exists $topk{$g}{$h}) {
                $v = 1; }
            elsif($g eq $h) { $v = 1; }
            print ODAT "\t$v"; }
        print ODAT "\n"; } }

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

