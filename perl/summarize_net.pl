#!/usr/bin/perl

# ===========================================================
# Name : summarize_net.pl
# Purpose : Provide #nodes, #edges and vertex-cover statitics
# Created : 27-09-2012
# Last Modified : Fri 28 Sep 2012 11:57:02 AM EDT
# Author(s) : Arjun Krishnan
# ===========================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, $idat, $idir, $ideg, $otab);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'idat=s' => \$idat,
            'idir' => \$idir,
          'ideg=i' => \$ideg,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $e);
my %all_edges = (); my %node_degree = (); my %node_edges = ();

print "\n$idat ...";
open DAT, "$idat" or die "Can't open $idat!";
while (<DAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    if($p[0] eq $p[1]) { next; } # Ignore self-edges
    if($p[2] == 0) { next; }     # Ignore edges with weight = 0

    # $node_degree{$p[0]}++; $node_degree{$p[1]}++;

    $e = join '__', sort($p[0], $p[1]);
    $all_edges{$e} = $p[2];
    $node_edges{$p[0]}{$e}++; $node_edges{$p[1]}{$e}++;
}
close DAT;

foreach my $n (keys %node_edges) {
    $node_degree{$n} = scalar keys %{$node_edges{$n}}; }

if($ideg) {
    # print "\n";
    foreach my $n (keys %node_degree) {
        if($node_degree{$n} < $ideg) { next; }

        foreach my $e (keys %{$node_edges{$n}}) {
            unless(exists $all_edges{$e}) { next; }
            delete $all_edges{$e}; }

        # print "\ndeleting $n ($node_degree{$n}) ...";

        delete $node_edges{$n};
        delete $node_degree{$n};
    }

    %node_edges = ();
    foreach my $e (keys %all_edges) {
        @p = split '__', $e;
        $node_edges{$p[0]}{$e}++; $node_edges{$p[1]}{$e}++; }

    %node_degree = ();
    foreach my $n (keys %node_edges) {
        $node_degree{$n} = scalar keys %{$node_edges{$n}}; }
}

my $num_nodes = scalar keys %node_degree;
my $num_edges = scalar keys %all_edges;

open TAB, ">$otab";
print TAB "No.Nodes\tFrac.Nodes\tNode.Degree\tNo.Edges\tFrac.Edges\tVC\n";

my %seen_edges = (); my $num_seene = my $num_seenn = my $vertex_cover = 0; 
my $frac_edges = my $frac_nodes = my $vc_measure = my $max_vcm = 0; my $print_max_vcm;
foreach my $n (sort {$node_degree{$b} <=> $node_degree{$a}} keys %node_degree) {
    $num_seenn++;
    foreach my $e (keys %{$node_edges{$n}}) {
        unless(exists $seen_edges{$e}) {
            $num_seene++;
            $seen_edges{$e}++; }
    }

    $frac_edges = ($num_seene / $num_edges);
    $frac_nodes = ($num_seenn / $num_nodes);
    $vc_measure = 2 * $frac_edges * (1 - $frac_nodes);
    $vc_measure /= $frac_edges + 1 - $frac_nodes;

    print TAB sprintf("%d\t%.5g\t%d\t%d\t%.5g\t%.5g\n", $num_seenn, $frac_nodes,
        $node_degree{$n}, $num_seene, $frac_edges, $vc_measure);

    if($vc_measure > $max_vcm) {
        $max_vcm = $vc_measure;
        $print_max_vcm = sprintf("%.5g", $vc_measure).' ( '.
            sprintf("%.5g", $frac_nodes).' nodes / '.sprintf("%.5g",
                $frac_edges).' edges )'; }

    if(($vertex_cover == 0) and ($num_seene == $num_edges)) {
        $vertex_cover = $num_seenn; }
}

print "\n\nNo. Nodes: $num_nodes\nNo. Edges: $num_edges";
print "\n\nVertex cover: $vertex_cover (", sprintf("%.5g",
    ($vertex_cover/$num_nodes));
print ") nodes\nmax VC_measure: $print_max_vcm\n\n";

__END__

=head1

Report summary statistics of input network.

=head1 USAGE

summarize_net.pl [--idat INPUT_DAT] [--otab OUTPUT_TAB] [--help]

=head1 DESCRIPTION

This script takes in a network in DAT format, assumes all edges with nonzero
wegihts to be valid edges, and then report back, i) no. nodes, ii) no. edges,
iii) greedy vertex cover, iv) vertex cover score, and v) cummulative dostribution of greedy vertex
cover.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input network in DAT format.

=item C<--idir>

(Optional) If provided, the network is tajen to be directed and two statistics
are provided -- for source and target nodes -- wherever appropriate. NOT
IMPLEMENTED YET.

=item C<--ideg>

(Optional) This option can be used after running the plain version and examining
the output table. If you notice that there are a few high-degree nodes that are
biasing the data, you can use this --ideg option to specify a maximum degree. If
provided, this cutoff will be used to retain only node (and hence the
respective) edges that have degree < ideg, and the summary will be recalculated.

=item C<--otab>

Output table.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Sep 27

=cut

