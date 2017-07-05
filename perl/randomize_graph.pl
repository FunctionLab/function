#!/usr/bin/perl

# ================================================
# Name : randomize_graph.pl
# Purpose : Randomize graph
# Created : 21-09-2012
# Last Modified : Tue 18 Dec 2012 05:58:43 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);

sub read_dat;
sub randomize_edges;
sub print_dat;

my ($help, $idat, $ineg, $odat);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'idat=s' => \$idat,
            'ineg' => \$ineg,
          'odat=s' => \$odat) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $time);

$time = runtime(); print "\n$time: Reading in $idat ...";
my $orig_edges = read_dat($idat);

$time = runtime(); print "\n$time: Randomizing graph ...";
my $perm_edges = randomize_edges($orig_edges);

$time = runtime(); print "\n$time: Printing $odat";
print_dat($perm_edges, $odat);

$time = runtime(); print "\n$time: DONE\n\n";

# Subroutines
# Read in a DAT
sub read_dat {
    my $dat = shift;
    my %edges = (); my (@q, $e);

    open DAT, "$dat" or die "Can't open $dat!";
    while (<DAT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @q = split '\t', $_;

        $e = join '__', sort($q[0], $q[1]);
        $edges{$e} = $q[2];
    }
    close DAT;

    return \%edges;
}

# Randomize graph
sub randomize_edges {
    my $edgeref = shift;
    my ($e1, $e2, @p1, @p2, $re1, $re2);

    my (%pos_edges, %rand_edges, @pose_array);
    if($ineg) {
        %pos_edges = %{$edgeref}; %rand_edges = ();
        @pose_array = shuffle keys %pos_edges; }
    else { %rand_edges = %{$edgeref}; }

    my $num_edges = scalar keys %{$edgeref};
    my $num_swaps = 0; my $tot_swaps = 3*$num_edges;
    $time = runtime(); print "\n\t$time: No. edges: $num_edges";

    if($ineg) {
        # ...
    }

    while($num_swaps <= $tot_swaps) {
        if($ineg) {
            $num_edges = scalar keys %pos_edges;
            $e1 = (keys %pos_edges)[int rand $num_edges];
            $e2 = (keys %pos_edges)[int rand $num_edges]; }
        else {
            $e1 = (keys %rand_edges)[int rand $num_edges];
            $e2 = (keys %rand_edges)[int rand $num_edges]; }

        if($e1 eq $e2) { next; }

        @p1 = split '__', $e1;
        @p2 = split '__', $e2;

        $re1 = $p1[0].'__'.$p2[1];
        $re2 = $p2[0].'__'.$p1[1];

        if($ineg) {
            if((exists $pos_edges{$re1}) or (exists $pos_edges{$re2})) { next; }
            if((exists $rand_edges{$re1}) or (exists $rand_edges{$re2})) { next; }
            $rand_edges{$re1} = $pos_edges{$e1}; delete $pos_edges{$e1};
            $rand_edges{$re2} = $pos_edges{$e2}; delete $pos_edges{$e2};
            $num_swaps++;
        }
        else {
            if((exists $rand_edges{$re1}) or (exists $rand_edges{$re2})) { next; }
            $rand_edges{$re1} = $rand_edges{$e1}; delete $rand_edges{$e1};
            $rand_edges{$re2} = $rand_edges{$e2}; delete $rand_edges{$e2};
            $num_swaps++;
        }

        unless($num_swaps % 1000) {
            $time = runtime(); print "\n\t$time: $num_swaps"; }
    }

    return \%rand_edges;
}

# Print edges into DAT
sub print_dat {
    my $edges = shift;
    my $dat = shift;

    my @q;

    open DAT, ">$dat";
    foreach my $e (keys %$edges) {
        @q = split '__', $e;
        print DAT "$q[0]\t$q[1]\t1\n";
    }
    close DAT;
}


__END__

=head1

Randomize input binary graph using edge-swap algorithm.

=head1 USAGE

randomize_graph.pl [--idat INPUT_DAT] [--odat OUTPUT_DAT] [--help]

=head1 DESCRIPTION

This script takes in a binary network in DAT format and randomizes the network
using the edge-swap algorithm to produce a permuted network preserving the
original node degree distribution.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input network in DAT format. Edges are taken to present between node pairs with
weights greater than 0, and absent between node pairs with weight 0 or not
present at all in the DAT file.

=item C<--ineg>

(Optional) If this option is provided, then during swapping, edges NOT in the
original network are created. This is useful for generating a negative standard
with the same degree distribution as that of the positives.

=item C<--odat>

Output network in DAT format. This will just be the list of edges with weights
equal to 1. Edges between all the other node pairs not present in this file are
to be deemed absent.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Sep 21

=cut

