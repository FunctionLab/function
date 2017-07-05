#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# my($imat1, $imat2, $omat, $idx1, $idx2) = @ARGV;
# if($#ARGV < 4) { $idx }
# elsif($#ARGV < 3) { $idx1 = $idx2 = 0; }
# imat1: Matrix 1
# imat2: Matrix 2
# omat: Combined matrix with columns of matrix 2 added to the right of columns of matrix1

my ($help, $imat1, $imat2, $ilc, $omat); my $idx1 = my $idx2 = 0; my $imaxg = 200;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<3);
GetOptions( 'help' => \$help,
         'imat1=s' => \$imat1,
         'imat2=s' => \$imat2,
          'idx1=i' => \$idx1,
          'idx2=i' => \$idx2,
             'ilc' => \$ilc,
         'imaxg=i' => \$imaxg,
          'omat=s' => \$omat) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open FH, "$imat1"; chomp(my @mat1=<FH>); close FH;
open GH, "$imat2"; chomp(my @mat2=<GH>); close GH;
open HH, ">$omat";

my %mat2=();
if($imat2 =~ /\.gmt$/) {
    my %gene_gs = ();
    foreach (@mat2) {
        my @p = split '\t', $_;
        my $gs = shift @p;
        (my $d = shift @p) =~ s/ \([0-9]*\)$//g;
        $d =~ s/[:;\/, ]/_/g; $d =~ s/[()'"]//g;
        $d =~ s/_-_/_/g; $d =~ s/__*/_/g;
        if((scalar @p) > $imaxg) { next; }
        foreach my $g (@p) {
            $gene_gs{$g}{$gs.':'.$d}++; } }

    foreach my $g (keys %gene_gs) {
        $mat2{$g} = join ' | ', keys %{$gene_gs{$g}}; } }
else {
    foreach (@mat2) {
        my @p = split '\t', $_;

        my $key2 = $p[$idx2]; if($ilc) { $key2 = lc($key2); }
        for(my $j=0; $j<=$#p; $j++) {
            if($j == $idx2) { next; }
            push(@{$mat2{$key2}}, $p[$j]); } } }

print HH "$mat1[0]\t"; my @col;
if($imat2 =~ /\.gmt$/) {
    (my $col = $imat2) =~ s/\.gmt$//g;
    $col =~ s/^.*\///g;
    print HH "$col\n"; }
else {
    @col = split '\t', $mat2[0]; splice(@col, $idx2, 1);
    print HH join "\t", @col; print HH "\n"; }

shift(@mat1);
foreach (@mat1) {
    my @p = split '\t', $_;

    print HH "$_";
    if($imat2 =~ /\.gmt$/) {
        unless(exists $mat2{$p[$idx1]}) { print HH "\t\n"; next; }
        print HH "\t$mat2{$p[$idx1]}\n"; }
    else {
        my $key1 = $p[$idx1]; if($ilc) { $key1 = lc($key1); }
        for(my $j=0; $j<=$#col; $j++) {
            print HH "\t";
            if(exists $mat2{$key1}) {
                print HH "${$mat2{$key1}}[$j]"; } }
        print HH "\n"; } }

    # if(exists $mat2{$p[$idx1]}) {
    # 	for(my $j=0; $j<=$#{$mat2{$p[$idx1]}}; $j++) {
    # 		print HH "\t${$mat2{$p[$idx1]}}[$j]"; } }
    # else { print HH "\t"; }

close HH;



__END__

=head1

Flexible code for binding and merging matrices by column keys.

=head1 USAGE

cbind_matrices.pl [--imat1 MATRIX1] [--idx1 COLUMN_WITH_KEYS_IN_MAT1] [--imat2
MATRIX2] [--idx2 COLUMN_WITH_KEYS_IN_MAT2] [--ilc] [--imaxg] [--omat
OUTPUT_MATRIX] [--help]

=head1 DESCRIPTION

This script takes in two matrices and combines them based on specified columns
in both inputs, adding the columns from mat2 along the corresponding matching
rows of mat1. The order of the matrices matters: mat1 merged with mat2 restricts
the output to only those matching keys/entries in mat1.

=head1 ARGUMENTS

=over 12

=item C<--imat1>

Input matrix 1 in tab-delimited format, with a header row containing column names.

=item C<--imat2>

Input matrix 2 in tab-delimited format, with a header row containing column names.

=item C<--omat>

Output file containing the bound/merged matrix.

=item C<--idx1>

(Optional) Index of column in matrix1 to use as reference for binding/merging.
Index begins at 0. Default is 0 (first column).

=item C<--idx2>

(Optional) Index of column in matrix2 to use as reference for binding/merging.
Index begins at 0. Default is 0.

=item C<--ilc>

(Optional) Converts keys in mat2 to lowercase before matching to keys in mat1.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

