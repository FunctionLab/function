#!/usr/bin/perl

# ================================================
# Name : combine_lists-into-matrix.pl
# Purpose : Combine all lists into a binary matrix
# Created : 20-12-2012
# Last Modified : Thu 01 Aug 2013 02:14:09 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my($help, @ilist, $iuni, $ival, $omat);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<2);
GetOptions( 'help' => \$help,
      'ilist=s{,}' => \@ilist,
          'iuni=s' => \$iuni,
            'ival' => \$ival,
          'omat=s' => \$omat) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if(scalar @ilist == 1) { die "Give me more than one list!\n"; }

my %all_genes = ();
open GH, "$iuni" or die "Can't open $iuni!";
while (<GH>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $all_genes{$_}++; }
close GH;

my @p;
open MAT, ">$omat"; print MAT "Gene";
my @list_array = (); my %gene_mem = (); my $l;
foreach my $list (@ilist) {
    $l = $list;
    $l =~ s/\.txt$//g; $l =~ s/^.*\///g; push(@list_array, $l);
    print MAT "\t$l";

    open LS, "$list" or die "Can't open $list!";
    while (<LS>) {
        if($_ =~ /^#/) { next; }
        chomp($_);
        if($ival) {
            if($_ =~ /\t/) { @p = split '\t', $_; }
            elsif($_ =~ / /) { @p = split ' ', $_; }
            unless(exists $all_genes{$p[0]}) { next; }
            $gene_mem{$l}{$p[0]} = $p[1]; }
        else {
            unless(exists $all_genes{$_}) { next; }
            $gene_mem{$l}{$_}++; } }
    close LS; }
print MAT "\n";

my $bit;
foreach my $g (sort keys %all_genes) {
    print MAT "$g";
    foreach my $l (@list_array) {
        $bit = 0; # if($ival) { $bit = 'NA'; }
        if(exists $gene_mem{$l}{$g}) {
            if($ival) {$bit = $gene_mem{$l}{$g}; }
            else { $bit = 1; } }
        print MAT "\t$bit"; }
    print MAT "\n"; }
close MAT;

__END__

=head1

Combine multiple lists (e.g. of genes) into a binary matrix of items x list-ids.

=head1 USAGE

./combine_lists-into-matrix.pl [--ilist INPUT_LISTS] [--o OUTPUT_MATRIX] [--help]

=head1 DESCRIPTION

This script takes in two or more lists of items (e.g. genes) and ouputs a matrix
of items along the rows and the list-identifiers along the columns, with 0 and 1
representing item membership in a list.

=head1 ARGUMENTS

=over 12

=item C<--ilist>

Input list files. Can provide many using wild-card matching.

=item C<--iuni>

List of union of all items in all the input lists.

=item C<--ival>

(Optional) If specified, the lists are assumed to contain one more column of
values for each entry, which will be recoreded and entered into the output
matrix.

=item C<--omat>

Output matrix file.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Dec 20

=cut

