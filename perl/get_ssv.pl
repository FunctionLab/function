#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);

my ($help, @idab, $istd, $inorm, $idck, $ossv); my $ifrac = 0.5;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 3);
GetOptions(  'help' => \$help,
        'idab=s{,}' => \@idab,
           'istd=s' => \$istd,
          'ifrac=f' => \$ifrac,
            'inorm' => \$inorm,
             'idck' => \$idck,
           'ossv=s' => \$ossv) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if($istd =~ /\.dab$/) { die "Gold-Std should be in DAT format."; }
my ($time, @p, $e, $c);

$time = runtime(); print "\n$time: Extracting data points ...\n";

my %dat_edges = (); my @dat_array = ();
my $maxf = 0; my ($dat, $d, $dck); my $ndat = 0;
foreach my $dab (@idab) {
    ($dat = $dab) =~ s/\.dab$/\.std\.dat/g; $dat =~ s/^.*\///g;
    `Dat2Dab -i $dab -e $istd -o $dat`;

    ($d = $dat) =~ s/\.std\.dat$//g; $d =~ s/^.*\///g;
    push(@dat_array, $d);

    $ndat++;
    $time = runtime(); print "\t$time:\t$ndat\t$dat ...\n";

    open DAT, "$dat" or die "Can't open $dat!";
    while (<DAT>) {
        @p = split '\t', $_;
        $e = join '__', sort($p[0], $p[1]);
        $dat_edges{$d}{$e} = $p[2];
        if($inorm and (abs($p[2]) > $maxf)) { $maxf = abs($p[2]); } }
    close DAT;

    `rm -f $dat`;

    if($idck) {
        ($dck = $istd) =~ s/\.dat$//g; $dck =~ s/^.*\///g;
        $dck = $d.'_'.$dck.'.dck';
        unless(-e $dck) {
            $time = runtime(); print "\t\t$time: dchecking ... $dck\n";
            `DChecker -w $istd -i $dab > $dck`; } } }


$time = runtime(); print "$time: Printing SSV ...\n";

my $npos = my $nneg = 0; my ($ne, $i);
open SSV, ">$ossv";
open STD, "$istd" or die "Can't open $istd!";
while (<STD>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    $e = join '__', sort($p[0], $p[1]);

    $ne = 0;
    foreach my $d (@dat_array) {
        if(exists $dat_edges{$d}{$e}) { $ne++; } }

    if(($ne/$ndat) < $ifrac) { next; }

    if($p[2] > 0) { $c = 1; $npos++; }
    else { $c = -1; $nneg++; }
    # $std_edges{$c}{$e}++;

    print SSV "$c"; $i = 1;
    foreach my $d (@dat_array) {
        if(exists $dat_edges{$d}{$e}) {
            if($inorm) {
                print SSV " $i:", sprintf("%.6g", ($dat_edges{$d}{$e}/$maxf)); }
            else { print SSV " $i:", sprintf("%.6g", $dat_edges{$d}{$e}); } }
        $i++; }
    print SSV "\n";
}
close STD;
close SSV;

print "\t+ve: $npos\t-ve: $nneg\n";
$time = runtime(); print "$time: DONE\n\n";


__END__

=head1

Prepares SSV for LIBLINEAR given gold-standard and feature datasets.

=head1 USAGE

get_ssv.pl [--istd GOLDSTD_DAT] [--idab FEATURES_DAT] [--ossv OUT_SSV] [--help]

=head1 DESCRIPTION

This script takes in a gold-standard of gene pairs and a bunch of feature
datasets, and assemble them into an SSV that LIBLINEAR can use.

=head1 ARGUMENTS

=over 12

=item C<--istd>

Gold standard gene-pairs in DAT format. Expected to contain edges with scores
0/1. If otherwise, those edges with scores > 0 will be deemed positive and the
rest will be deemed negative.

=item C<--idab>

(Optional) Input feature datasets in DAB format. Multiple datasets can be provided using
wildcard matching. If not provided (by default), all DABs in the current working
directory will be used. Each DAB is expected to be a complete matrix with
non-edges. If not, use Dat2Dab with -Z option to create a complete matrix.

=item C<--ifrac>

(Optional) Fraction of datasets a gold-std edge has to have values in to be
included in the SSV. Default 0.5.

=item C<--inorm>

(Optional) Normalize feature values to be b/w 0 and 1.

=item C<--idck>

(Optional) Evaluate each input dataset against the gold-std using DChecker.

=item C<--ossv>

Output file in SSV format.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 Jan 29

=cut

