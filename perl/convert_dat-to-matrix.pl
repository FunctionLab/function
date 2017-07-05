#!/usr/bin/perl

# ================================================
# Name : convert_dat-to-matrix.pl
# Purpose : Convert dat to a gene-gene matrix
# Created : 14-08-2013
# Last Modified : Tue 30 Dec 2014 03:34:15 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, @idat, $irowg, $icolg, $icut);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(  'help' => \$help,
        'idat=s{,}' => \@idat,
          'irowg=s' => \$irowg,
          'icolg=s' => \$icolg,
           'icut=f' => \$icut,   ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my %rowg = my %colg = ();

if($irowg) {
    open RF, "$irowg" or die "Can't open $irowg!";
    while (<RF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $rowg{$_}++; }
    close RF; }

if($icolg) {
    open CF, "$icolg" or die "Can't open $icolg!";
    while (<CF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $colg{$_}++; }
    close CF; }

foreach my $dat (@idat) {
    print "\n$dat";
    my %val = my %hrow = my %hcol = ();
    my %gene_ne = my %gene_mean = ();

    open DAT, "$dat" or die "Can't open $dat!";
    while (<DAT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        
        my $v = $p[2];
        if($icut) { if($v >= $icut) { $v = 1; } else { $v = 0; } }

        $gene_ne{$p[0]}++;
        unless(exists $gene_mean{$p[0]}) { $gene_mean{$p[0]} = 0; }
        $gene_mean{$p[0]} += ($v - $gene_mean{$p[0]})/$gene_ne{$p[0]};

        $gene_ne{$p[1]}++;
        unless(exists $gene_mean{$p[1]}) { $gene_mean{$p[1]} = 0; }
        $gene_mean{$p[1]} += ($v - $gene_mean{$p[1]})/$gene_ne{$p[1]};
        
        my $e = join '__', sort($p[0], $p[1]);
        $val{$e} = $v;

        if($irowg) {
            if(exists $rowg{$p[0]}) { $hrow{$p[0]}++; }
            if(exists $rowg{$p[1]}) { $hrow{$p[1]}++; } }
        else {
            $hrow{$p[0]}++; $hrow{$p[1]}++; }
        if($icolg) {
            if(exists $colg{$p[0]}) { $hcol{$p[0]}++; }
            if(exists $colg{$p[1]}) { $hcol{$p[1]}++; } }
        else {
            $hcol{$p[0]}++; $hcol{$p[1]}++; } }
    close DAT;

    my @arow = sort keys %hrow;
    my @acol = sort keys %hcol;

    (my $mat = $dat) =~ s/^.*\///g; $mat =~ s/\.dat$/\.mat/g;
    open MAT, ">$mat";

    print MAT "Gene";
    for(my $j=0; $j<=$#acol; $j++) {
        print MAT "\t$acol[$j]"; }
    print MAT "\n";

    for(my $i=0; $i<=$#arow; $i++) {
        print MAT "$arow[$i]";
        for(my $j=0; $j<=$#acol; $j++) {
            if($arow[$i] eq $acol[$j]) {
                printf MAT "\t%.6g", $gene_mean{$arow[$i]}; }
            else {
                my $e = join '__', sort($arow[$i], $acol[$j]);
                print MAT "\t$val{$e}"; } }
        print MAT "\n"; }
    close MAT; }
print "\n\n";


__END__

=head1

Convert one or more DAT files into respective gene-gene association matrices.

=head1 USAGE

convert_dat-to-mat.pl [--idat INPUT_DAT] [--irowg ROW-GENE_FILTER] [--icolg
COL-GENE_FILTER] [--icut SCORE_CUTOFF] [--help]

=head1 DESCRIPTION

This script takes in one or more networks in DAT format and converts each into a
gene association matrix. The genes that occupy the rows and columns can be
filtered.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input DAT files that need to be converted to matrices. More than one can be
supplied usig wild-cards.

=item C<--irowg>

(Optional) List of genes to retain as rows.

=item C<--icolg>

(Optional) List of genes to retain as columns.

=item C<--icut>

(Optional) Edge weight cutoff to use. Edges above this cutoff will be recorded
as 1 and the other as 0.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

