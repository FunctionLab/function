#!/usr/bin/perl

# ================================================
# Name : sample_weighted-gold-std.pl
# Purpose : Get a random sample of edges from equally spaced weight bins.
# Created : 15-11-2012
# Last Modified : Thu 15 Nov 2012 02:01:18 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my ($help, $igs, $ogs); my $inbin = 10;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions( 'help' => \$help,
           'igs=s' => \$igs,
         'inbin=i' => \$inbin,
           'ogs=s' => \$ogs) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my ($time, @p);
print "\nFiltering edges in bins ... ";

my $min_edges = 10000; my ($str, $end, $odat, $nume); my %subsets = ();
for(my $i=0; $i<0.999; $i+=(1/$inbin)) {
    $str = sprintf("%.3f", $i);
    $end = sprintf("%.3f", $i+(1/$inbin)-0.001);
    if($end == 0.999) { $end = sprintf("%.3f", 1);  }
    print "\n\t$str - $end";
    
    $odat = 'gs_'.$str.'-'.$end.'.dat';
    $subsets{$str.'-'.$end} = $odat;

    `Filterer -i $igs -o $odat i$str=$end`;
    `perl -MList::Util -e 'print List::Util::shuffle <>' $odat > temp.dat`;
    `mv temp.dat $odat`;
    
    chomp($nume = `wc -l $odat`); $nume =~ s/ $odat$//g;
    if($nume < $min_edges) { $min_edges = $nume; }

    $time = runtime(); print "\t$time\t$odat - $nume";
}

print "\nMin edges: $min_edges\nGetting subsets ...";

`rm -f $ogs`;
foreach (sort keys %subsets) {
    $time = runtime(); print "\n\t$time: $_";
    `head -$min_edges $subsets{$_} >> $ogs`; }

# $ogs =~ s/\.dat/\.over-sampled\.dat/g;

# open OGS, "$ogs";
# while (<OGS>) {
#     chomp($_); @p = split '\t', $_;
#     
#     for(my $r=0; $r<=$#f; $i++) {
#     }
# }
# close OGS;

print "\nDONE\n\n";

__END__

=head1

Sample weighted gold-std edges with equal numbers across bins.

=head1 USAGE

sample_weighted-gold-std.pl [--igs INPUT_GS_DAT/DAB] [--inbin NUM_BINS] [--ogs
OUTPUT_GS_DAT] [--help]

=head1 DESCRIPTION

This script takes in a collection of gold-std edges with weights (between 0 and
1) and samples equal number of edges from uniformly-spaced bins. The number of
bins can be specified by the user, and the minimum number of edges within on of
the bins is used as the number of edges to sample across all bins.

=head1 ARGUMENTS

=over 12

=item C<--igs>

Input weighted-gold-std file in DAT/DAB format.

=item C<--inbin>

No. of bins to split the range between 0 and 1.

=item C<--ogs>

Output sampled gold-std file in DAT format.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Nov 15

=cut

