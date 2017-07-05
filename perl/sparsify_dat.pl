#!/usr/bin/perl

# ================================================
# Name : sparsify_dat.pl
# Purpose : Sparsify network by applying edge-cutoff and retaining only top k neighbors of genes
# Created : 06-01-2015
# Last Modified : Mon 20 Feb 2017 08:05:49 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);


my ($help, $idat, $icut, $ikeeps, $odat); my $itopk = 50;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions(     'help' => \$help,
              'idat=s' => \$idat,
             'itopk=i' => \$itopk,
              'ikeeps' => \$ikeeps,
              'icut=f' => \$icut,
              'odat=s' => \$odat    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my $time = runtime(); print "\n$time: Parsing DAT ...";

my %net = (); my $one = 0;
open IDAT, "$idat" or die "Can't open $idat!";
while (<IDAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $one++;

    $net{$p[0]}{$p[1]} = $p[2];
    $net{$p[1]}{$p[0]} = $p[2]; }
close IDAT;

my @agenes = sort keys %net;
my $ong = scalar @agenes;


$time = runtime(); print "\n$time: Indexing top k genes per gene ...";

my %edges = my %all_topk = my %topk = ();
foreach my $g (@agenes) {
    my @a = sort {$net{$g}{$b} <=> $net{$g}{$a}} keys %{$net{$g}};

    my $k = $itopk;
    if((scalar @a) < $itopk) { $k = scalar @a; }

    foreach my $h (@a[0..($k-1)]) {
        if($icut and ($net{$g}{$h} < $icut)) {
            next; }

        $all_topk{$h}++;
        $topk{$g}{$h}++;

        my $e = join '__', sort($g, $h);
        $edges{$e} = $net{$g}{$h}; } }


my @tgenes = sort keys %all_topk;

$time = runtime(); print "\n$time: Printing output file ...";

open ODAT, ">$odat";
if($odat =~ /\.pcl$/) {
    foreach (@tgenes) {
        print ODAT "\t$_"; }
    print ODAT "\n";

    foreach my $g (@agenes) {
        print ODAT "$g";

        my $ntopk = 0;
        foreach my $h (@tgenes) {
            my $v = 0;
            if(exists $topk{$g}{$h}) {
                $v = 1; $ntopk++; }
            #elsif($g eq $h) { $v = 1; }
            print ODAT "\t$v"; }
        print ODAT "\n";

        print "$g\t$ntopk\n"; }
    close ODAT;

    $odat =~ s/\.pcl$/\.dat/g; }


open ODAT, ">$odat";
my $sne = 0; my %genes = ();
foreach my $e (keys %edges) {
    my @p = split '__', $e;
    if((exists $topk{$p[0]}{$p[1]}) and (exists $topk{$p[1]}{$p[0]})) {
        $sne++;
        $genes{$p[0]}++; $genes{$p[1]}++;

        if($odat =~ /\.dat$/) {
            print ODAT "$p[0]\t$p[1]\t$edges{$e}\n"; } } }
close ODAT;

my $sng = scalar keys %genes;


print "\n\n$idat\n\tOrig: $ong genes\t$one edges\n\tSpar: $sng genes\t$sne edges";
print "\t", scalar @tgenes, " top genes\n\n";


__END__

=head1

Sparsify network to retain top edges above cutoff.

=head1 USAGE

sparsify_dat.pl [--idat NETWORK_DAT] [--itopk TOP-K_NEIGHBORS] [--icut
EDGE_CUTOFF] [--odat NETWORK_DAT] [--help]

=head1 DESCRIPTION

This script takes in a network in DAT format, and sparsifies it by retaining
only the top k edges per node, and further removing edges below a specified
cutoff (if provided).

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input network in DAT format.

=item C<--itopk>

(Optional) Number of top edges to retain per node. Default is 50.

=item C<--icut>

(Optional) Edge weight cutoff.

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

