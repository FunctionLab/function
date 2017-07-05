#!/usr/bin/perl

# ================================================
# Name : calc_network-geneset-connectivity_part1.pl
# Purpose : Calculate connectivity of real and random genesets within network
# Created : 24-01-2013
# Last Modified : Sun 11 Jan 2015 07:30:12 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);


my $hubber = '/Genomics/ogtr04/arjunk/bin/Hubber';


my ($help, $igmt, @iglist, $idab);
my $iming = 5; my $imaxg = 100;
my $iperm = 5000; my $inetbg;

pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 3);
GetOptions( 'help' => \$help,
          'idab=s' => \$idab,
          'igmt=s' => \$igmt,
          'inetbg' => \$inetbg,
     'iglist=s{,}' => \@iglist,
         'iming=i' => \$iming,
         'imaxg=i' => \$imaxg,
         'iperm=i' => \$iperm    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


# Index network genes
(my $ideg = $idab) =~ s/\.dab$/\.deg/g; $ideg =~ s/^.*\///g;
unless(-e $ideg) {
    `Dat2Dab -i $idab -P | tail -n +3 > $ideg`; }

my %netg = ();
open GEN, "$ideg" or die "Can't open $ideg!";
while (<GEN>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $netg{$p[0]} = $p[1]; }
close GEN;

my @all_netg = sort keys %netg;


# Index geneset collection
my ($gs_genes, $gs_size, $all_gsg, $N);

if($igmt) {
    undef @iglist;
    ($gs_genes, $gs_size, $all_gsg) = parse_gmt($igmt, \%netg);
    $N = scalar @$all_gsg; }
# else { }  ... # use this to index individual gene lists

if($inetbg) {
    $N = @all_netg; }


# Index genesets by size and generate gene lists
(my $gdir = $igmt) =~ s/\.gmt$/_gen/g; $gdir =~ s/^.*\///g;
(my $hdir = $igmt) =~ s/\.gmt$/_hub/g; $hdir =~ s/^.*\///g;
`mkdir -p $gdir`; `mkdir -p $hdir`;

my %size_gs = ();
foreach my $gs (sort keys %$gs_size) {
    if(($gs_size->{$gs} < $iming)
            or ($gs_size->{$gs} > $imaxg)) { next; }

    (my $gsid = $gs) =~ s/:/__/g;
    print "\n$gsid";

    my $gfile = $gdir.'/'.$gsid;
    unless(-e $gfile) {
        open GS, ">$gfile";
        foreach my $g (sort keys %{$gs_genes->{$gs}}) {
            print GS "$g\n"; }
        close GS; }

    my $hfile = $hdir.'/'.$gsid;
    unless(-e $hfile) {
        `qsub -m n -N $gsid.hub -l 1hr -cwd "$hubber -i $idab $gfile > $hfile"`; }

    $size_gs{$gs_size->{$gs}}{$gs}++; }


# Printing random genesets
my $jcount;
foreach my $sz (sort keys %size_gs) {
    for(my $i=1; $i<$iperm; $i++) {
        chomp($jcount = `qstat -u arjunk | wc -l`);

        while($jcount > 500) {
            print "\n\t$jcount: sleeping for a min...";
            sleep 60;
            chomp($jcount = `qstat -u arjunk | wc -l`); }

        my $gsid = 's'.$sz.'.r'.sprintf("%05d", $i);
        print "\n$gsid";

        my $gfile = $gdir.'/'.$gsid;
        unless(-e $gfile) {
            my @rand = shuffle 0..($N-1);
            if($inetbg) {
                @rand = @all_netg[@rand[0..($sz-1)]]; }
            else {
                @rand = @$all_gsg[@rand[0..($sz-1)]]; }

            open GS, ">$gdir/$gsid";
            foreach my $g (@rand) {
                print GS "$g\n"; }
            close GS; }

        my $hfile = $hdir.'/'.$gsid;
        unless(-e $hfile) {
            `qsub -m n -N $gsid.hub -l 1hr -cwd "$hubber -i $idab $gfile > $hfile"`; } } }
print "\n\n";



# Subroutines
sub parse_gmt {
    my $gmt = shift;
    my $bg = shift;

    my %genes = my %size = my %allg = ();

    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @q = split '\t', $_;

        my $gs = shift @q;

        $q[0] =~ s/ \([0-9]*\)$//g;
        shift @q;

        foreach my $g (@q) {
            unless(exists $bg->{$g}) { next; }
            $genes{$gs}{$g}++;
            $allg{$g}++;
            $size{$gs}++; } }
    close GMT;

    my @allg = sort keys %allg;

    return(\%genes, \%size, \@allg); }


__END__

=head1

Calculate connectivity of genesets within network and evaluate significance.

=head1 USAGE

./gsea_gmt_NPAGE.pl [--igmt GENESETS_GMT] [--iglist GLISTS_TXT] [--idat NETWORK_DAT] [--o OUT_TAB] [--help]

=head1 DESCRIPTION

This script takes in a collection of genesets and calculates their connectivity
in a network provided by the user. The siginficance of the connectivity is
tested using a permutation test.

=head1 ARGUMENTS

=over 12

=item C<--igmt>

Input collection of genesets in GMT format. Use this or the --iglist option.

=item C<--iglist>

Input collection of genesets as a collection of text files matched using
wild-card. Each file should list genes one-per-line. Provide either --igmt or
--iglist. If --igmt is provided, --iglist is ignored.

=item C<--idab>

Input network in DAT/DAB format.

=item C<--otab>

Output table.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 Jan 24

=cut

