#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);
use Graph::Undirected;

my ($help, $in_snet, $in_qnet, $in_g1, $in_g2, $in_pspl, $out_net, $time);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
         'isnet=s' => \$in_snet,
         'iqnet=s' => \$in_qnet,
           'ig1=s' => \$in_g1,
           'ig2=s' => \$in_g2,
            'ispl' => \$in_pspl,
            'onet' => \$out_net) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

$time = runtime(); print "$time: Indexing subject network ...\n";

open SNET, "$in_snet" or die "Can't open $in_snet!";
my $snet = Graph::Undirected->new;
my %snet_edges = (); my %snet_degree = (); my %snet_neighbors = ();
my (@p, $edge);
while (<SNET>) {
    chomp($_);
    @p = split '\t', $_;
    $edge = join '__', sort($p[0], $p[1]);
    if((exists $snet_edges{$edge}) or ($p[2] <= 0)) { next; }

    $snet_edges{$edge} = $p[2];
    $snet_degree{$p[0]} += $p[2]; $snet_degree{$p[1]} += $p[2];
    $snet_neighbors{$p[0]}{$p[1]}++; $snet_neighbors{$p[1]}{$p[0]}++;
    $snet->add_weighted_edge($p[0], $p[1], $p[2]);
}
close SNET;

my $tot_genes = scalar(keys %snet_degree);

$time = runtime(); print "$time: Calculating node and edge weights ...\n";

my $wsnet = Graph::Undirected->new;
my %wsnet_degree = ();
my ($w_uv, $p_u, $p_v); # my $tot_wedges = 0;
foreach my $edge (keys %snet_edges) {
    @p = split '__', $edge;
    $p_u = ($snet_degree{$p[0]}/$tot_genes);
    $p_v = ($snet_degree{$p[1]}/$tot_genes);
    
    if($in_pspl) { $w_uv = sqrt($p_u*$p_v); }
    else { $w_uv = sqrt((1-$p_u)*(1-$p_v)); }

    $wsnet_degree{$p[0]} += $w_uv; $wsnet_degree{$p[1]} += $w_uv;
    $wsnet->add_weighted_edge($p[0], $p[1], $w_uv);
}

if($in_g1) {
    open GS, "$in_g1" or die "Can't open $in_g1!";
    my %gs1 = (); while(<GS>) { $gs1{chomp($_)}++; } close GS;
}

if($in_g2) {
    open GS, "$in_g2" or die "Can't open $in_g2!";
    my %gs2 = (); while(<GS>) { $gs2{chomp($_)}++; } close GS;
}

if($in_qnet) {
}


__END__

=head1

Calculation of topological-overlap-measure (TOP; Default) or shortest-path-length (SPL)
based on an underlying network.

=head1 USAGE

./get_tom_spl.pl [--isnet SBJCT_NETWORK_DAT] [--iqnet QUERY_NETWORK_DAT] [--ig1
GENE_LIST1] [--ig2 GENE_LIST2] [--ispl] [--onet OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in a subject network and i) query network, or ii) geneset(s), and
calculates TOM (default) or SPL for i) all the edges in the query network, or ii) all
pairs of genes within a geneset or spanning the two genesets. On the otherhand,
if only the subject network is provided (without the other file input options),
then, based on the chosen paramter, the TOM or SPL measures are calculated
for all pairs of genes in that network.

In all these calculations, for a pair of genes under consideration at any point,
the existence of a direct edge between them in the subject network is actively
ignored. Therefore, the calculations are a measure of either how many 'true'
neghbors there are between the genes (TOM), or in how many steps can you reach
one gene from the other in the absence of a direct edge.

=head1 ARGUMENTS

=over 12

=item C<--isnet>

Subject network in DAT format (Required). Only non-zero
edges will do. Edges with zero or neagtive weights are considered to be absent
edges.

=item C<--iqnet>

Query network in DAT format (Optional).

=item C<--ig1>

Gene list 1 (Optional). If provided without 'ig2', all pairs of genes in this
geneset are considered.

=item C<--ig2>

Gene list 2 (Optional; but can be provided only along with 'ig1'). Gene pairs
spanning the 'ig1' & 'ig2' genesets are considered.

=item C<--ispl>

Calcualte SPL instead of TOM (Optional).

=item C<--onet>

Output network file in DAT format (Required).

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 24

=cut

