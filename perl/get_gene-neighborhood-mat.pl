#!/usr/bin/perl

# ================================================
# Name : get_gene-neighborhood-mat.pl
# Purpose : Get a matrix of input genes along columns and network genes along rows
# Created : 01-10-2014
# Last Modified : Wed 01 Oct 2014 01:36:52 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $dat2dab = '/r04/arjunk/bin/Dat2Dab';


my ($ihelp, $idab, $igenes, $icut, $omat);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'idab=s' => \$idab,
        'igenes=s' => \$igenes,
          'icut=f' => \$icut,
          'omat=s' => \$omat    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


(my $inetg = $idab) =~ s/\.dab$/\.genes/g;
`$dat2dab -i $idab -E > $inetg`;

my %netg = ();
open NETG, "$inetg" or die "Can't open $inetg!";
while (<NETG>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $netg{$_}++; }
close NETG;


my %genes = ();
open GEN, "$igenes" or die "Can't open $igenes!";
while (<GEN>) {
    if($_ =~ /^#/) { next; }
    chomp($_); unless(exists $netg{$_}) { next; }
    $genes{$_}++; }
close GEN;



__END__

=head1

Generates a gene-neighborhood matrix with input genes along columns and network
genes along rows.

=head1 USAGE

get_gene-neighborhood-mat.pl [--idab NETWORK_DAT/DAB] [--igenes LIST_OF_GENES]
[--icut EDGE_CUTOFF] [--omat OUTPUT_MAT] [--help]

=head1 DESCRIPTION

This script takes a network in DAT/DAB format and a list of genes (one-per-line)
and generates a matrix with input genes along the columns, all network genes
along the rows and the correspodning network edge scores in the cells of the
matrix.

=head1 ARGUMENTS

=over 12

=item C<--idab>

Input network in DAT/DAB format.

=item C<--igenes>

Input list of genes.

=item C<--icut>

(Optional) Cutoff on the edges. Edges with score < cutoff will be set to 0.

=item C<--omat>

Output matrix of gene neighborhoods.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2014 Oct 01

=cut

