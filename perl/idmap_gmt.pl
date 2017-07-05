#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);


my ($help, $igmt, $imap, $ifiltgs, $ifiltg, $itrimr, $ilogm, $ogmt);
my $imeasr = 'median';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(     'help' => \$help,
              'igmt=s' => \$igmt,
              'imap=s' => \$imap,
              'itrimr' => \$itrimr,
            'imeasr=s' => \$imeasr,
               'ilogm' => \$ilogm,
            'ifiltg=s' => \$ifiltg,
            'ifiltgs=s' => \$ifiltgs,
              'ogmt=s' => \$ogmt        ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %oid2nid = ();
open MAP, "$imap" or die "Can't open $imap!";
while (<MAP>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    if(($p[0] eq '') or ($p[0] =~ /\?/) or ($p[0] eq 'NA')) { next; }
    if(($p[1] eq '') or ($p[1] =~ /\?/) or ($p[1] eq 'NA')) { next; }

    my $oid = shift @p;
    if($itrimr) { $oid =~ s/\.[0-9][0-9]*$//g; }

    foreach my $nid (@p) {
        $oid2nid{$oid}{$nid}++; } }
close MAP;


my %filtgs = ();
if($ifiltgs) {
    open FGS, "$ifiltgs" or die "Can't open $ifiltgs!";
    while (<FGS>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $filtgs{$p[0]}++; }
    close FGS; }


my %filtg = ();
if($ifiltg) {
    open FG, "$ifiltg" or die "Can't open $ifiltg!";
    while (<FG>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $filtg{$p[0]}++; }
    close FG; }


open OGMT, ">$ogmt";
open IGMT, "$igmt" or die "Can't open $igmt!";
while (<IGMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $gs = shift @p;
    (my $desc = shift @p) =~ s/ \([0-9]*\)$//g;

    my %genes = ();
    foreach my $oid (@p) {
        foreach my $nid (keys %{$oid2nid{$oid}}) {
            $genes{$nid}++; } }

    print OGMT "$gs\t$desc (", (scalar keys %genes), ")";
    foreach my $nid (keys %genes) {
        print OGMT "\t$nid"; }
    print OGMT "\n"; }
close IGMT;
close OGMT;


# Calcualte mac
sub max {
    my $aref = shift;

    my @array = sort {$b <=> $a} grep {!/^NA$/} @$aref;
    return $array[0]; }


# Calculate mean
sub mean {
    my $aref = shift;
    
    if((scalar @$aref) == 1) {
        return ${$aref}[0]; }

    else {
        my $mean = my $n = 0;
        foreach (@$aref) {
            $n++;
            my $v = $_; if($ilogm) { $v = log($v + 0.0001); }
            $mean += ($v - $mean)/$n; }

        if($ilogm) { return (exp($mean) - 0.0001); }
        else { return $mean; } } }


# Calculates median
sub median {
    my $aref = shift;

    if((scalar @$aref) == 1) {
        return ${$aref}[0]; }

    else {
        my @array = sort {$a <=> $b} grep {!/^NA$/} @$aref;
        if(@array % 2) { return $array[@array/2]; }
        else { return (($array[(@array/2)-1] + $array[@array/2]) / 2); } } }


# Calculate sum
sub sum {
    my $aref = shift;

    my $sum = 0;
    foreach (@$aref) { $sum += $_; }

    return $sum; }


__END__

=head1

Map simple matrix from one row-identifier space to another.

=head1 USAGE

idmap_mat.pl [--igmt INPUT_MATRIX] [--imap ID-MAP] [--ifiltg ROWS-TO-FILTER]
[--ifiltgs COLUMNS-TO-FILTER] [--ogmt OUTPUT_MATRIX] [--help]

=head1 DESCRIPTION

This script will map a matrix from one row id-space to another. When there is
many-to-many mapping, all the old ids that map to multiple new ids will be
removed, and when multiple given ids map to a single new id, the rows of the old
ids will be averaged column-wise and assigned to the new id.

=head1 ARGUMENTS

=over 12

=item C<--igmt>

Input matrix in a simple tab-separated format with header row containing column
identifiers, first columns containing row identifiers, and the rest of the file
filled with values.

=item C<--imap>

Mapping of row identifiers from one idspace to another. Tab-separated values of
<id1> <ide2>. Could contain many-to-many mapping, in one or many lines.

=item C<--ifiltg>

(Optional) List of row identifiers (in the new id space) that should be retained
in the output file.

=item C<--ifiltgs>

(Optional) List of column identifiers that should be retained in the output file.

=item C<--itrimr>

(Optional) Trim row identifiers off numerical suffixes before mapping
(e.g. ENSG00234234.1 to ENSG00234234).

=item C<--ilogm>

(Optional) Calculate mean in log-space and unlog in the output.

=item C<--ogmt>

Output matrix in the new row-identifier space.

=item C<--help>

Prints this help message.

=back

=head1 PARAMETERS

=over 12

=item C<--p>

Parameter

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

