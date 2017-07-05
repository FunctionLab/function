#!/usr/bin/perl

# ================================================
# Name : get_gene-neighbor-genesets.pl
# Purpose : Get genesets from neighbors of genes in network
# Created : 26-09-2013
# Last Modified : Thu 26 Sep 2013 10:12:13 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);


my (@idat, $help, $time); my $icut = 1; my $inbh = 1;
my $isym = '/home/arjunk/data/mappings/human_genes_id-symbol-description.txt';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
       'idat=s{,}' => \@idat,
          'inbh=i' => \$inbh,
          'icut=f' => \$icut ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %gene_sym = ();
open SYM, "$isym" or die "Can't open $isym!";
while (<SYM>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $gene_sym{$p[0]} = $p[1]; }
close SYM;


foreach my $dat (@idat) {
    $time = runtime(); print "\n$time: $dat";

    my %gene_neighbors = ();
    open DAT, "$dat" or die "Can't open $dat!";
    while (<DAT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        if($p[2] < $icut) { next; }

        $gene_neighbors{$p[0]}{$p[1]}++;
        $gene_neighbors{$p[1]}{$p[0]}++; }
    close DAT;

    (my $gmt = $dat) =~ s/\.dat/\.gmt/g;
    open GMT, ">$gmt";
    foreach my $g (keys %gene_neighbors) {
        my $nnbh = scalar keys %{$gene_neighbors{$g}};
        if($nnbh < $inbh) { next; }
        print GMT "$g\t";
        if(exists $gene_sym{$g}) { print GMT "$gene_sym{$g}"; }
        else { print GMT "$g"; }
        print GMT " ($nnbh)";

        foreach my $ng (keys %{$gene_neighbors{$g}}) {
            print GMT "\t$ng"; }
        print GMT "\n"; }
    close GMT; }

$time = runtime(); print "\n\n$time: DONE\n\n";


__END__

=head1

Generates a collection of neighborhood-genesets for each gene in a network.

=head1 USAGE

get_gene-neighbor-genesets.pl [--idat NETWORK_DAT] [--icut VALUE_DOUBLE] [--inbh NUM_NEIGHBORS_INT] [--help]

=head1 DESCRIPTION

This script takes in one or more network files in DAT format, and for each,
generates a collection of genestes in GMT format where each geneset corresponds
to a gene and the set is all the network neighbors of that gene passing a
threshold supplied through the --icut parameter.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input network file in DAT format. More than one can be supplied using a
wild-cards.

=item C<--icut>

(Optional) Only edge score >= this threshold will be used to select neighbors of
genes. Default is 1.

=item C<--inbh>

(Optional) No. of neighbors a gene needs to have to be printed in the GMT.
Default is 1.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 Sep 26

=cut

