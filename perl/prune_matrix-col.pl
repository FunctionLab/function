#!/usr/bin/perl

# ================================================
# Name : prune_matrix-col.pl
# Purpose : 
# Created : 20-12-2012
# Last Modified : Fri 21 Dec 2012 02:50:04 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my($help, $imat, $iclu, $omat);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'iclu=s' => \$iclu,
          'omat=s' => \$omat) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my @p;

my %cluster_col = (); my %uniq_col = ();
open CLU, "$iclu" or die "Can't open $iclu!";
while (<CLU>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    foreach my $c (@p) { $cluster_col{$c}++; }
    $uniq_col{$p[int(rand(scalar @p))]}++; } # choosing one per randomly per cluster
close CLU;

open IMAT, "$imat" or die "Can't open $imat!";
chomp(my @mat = <IMAT>); close IMAT;
my @all_col = split '\t', shift @mat;

my $tot_col = my $nclu_col = my $nuniq_col = 0;
open OMAT, ">$omat"; print OMAT "$all_col[0]";
for(my $j=1; $j<=$#all_col; $j++) {
# foreach (@all_col) {
    $tot_col++;
    if(exists $cluster_col{$all_col[$j]}) {
        $nclu_col++;
        if(exists $uniq_col{$all_col[$j]}) {
            $nuniq_col++;
            print OMAT "\t$all_col[$j]"; } }
    else { print OMAT "\t$all_col[$j]"; }
}
print OMAT "\n";

print "\nTot. col: $tot_col\nNo. clustered col: $nclu_col\nNo. unique col: $nuniq_col\n\n";

foreach (@mat) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    
    print OMAT "$p[0]";
    for(my $j=1; $j<=$#p; $j++) {
        if(exists $cluster_col{$all_col[$j]}) {
            if(exists $uniq_col{$all_col[$j]}) {
                print OMAT "\t$p[$j]"; } }
        else { print OMAT "\t$p[$j]"; }
    }
    print OMAT "\n";
}
close OMAT;

__END__

=head1

Prune matrix columns to get one column per cluster of similar columns.

=head1 USAGE

./prune_matrix-col.pl [--imat INPUT_MATRIX] [--p INPUT_COL-CLUSTERS] [--o
OUTPUT_MATRIX] [--help]

=head1 DESCRIPTION

This script takes in an input matrix and a clustering output of very-similar
columns, and outputs a pruned matrix with only one column per cluster of similar
columns.

=head1 ARGUMENTS

=over 12

=item C<--imat>

Input matrix.

=item C<--iclu>

SPICi cluster output.

=item C<--omat>

Output matrix.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

