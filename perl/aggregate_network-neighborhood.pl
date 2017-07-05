#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);

my ($help, $iglist, $ibg, @inet, $time, @p);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
        'iglist=s' => \$iglist,
           'ibg=s' => \$ibg,
       'inet=s{,}' => \@inet) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my %qgenes = ();
open GH, "$iglist" or die "Can't open $iglist!";
while(<GH>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $qgenes{$_}++;
}
close GH;

print "\nExtracting network connections ...";
my ($odat, $gene2); my %gene_val = (); my @nets = ();
foreach my $net (@inet) {
    # if($net =~ /global\.dab/) { next; }
    # if($net =~ /tissues_cell_types/) { next; }
    push(@nets, $net); $nets[$#nets] =~ s/^.*\///g; $nets[$#nets] =~ s/\.q*dab//g;
    print "\n$nets[$#nets] ...";

    foreach my $gene (keys %qgenes) {
        ($odat = $net) =~ s/^.*\///g; $odat =~ s/\.q*dab/\.dat/g;
        $odat = $gene.'-'.$odat;
        $time = runtime(); print "\n\t$time:\t$gene";

        `Dat2Dab -i $net -l $gene > $odat`;

        open DAT, "$odat"; $odat =~ s/^$gene-//g; $odat =~ s/\.dat$//g;
        while (<DAT>) {
            chomp($_); @p = split '\t', $_;
            $gene2 = $p[1]; if($p[1] eq $gene) { $gene2 = $p[0]; }
            push(@{$gene_val{$gene}{$gene2}}, $p[2]);
        }
        close DAT;

        `rm -f $odat`;
    }
}

my ($omat, @other_genes); print "\n\nPrinting matrices ...";
foreach my $gene (keys %qgenes) {
    $time = runtime(); print "\n\t$time: $gene ...";

    $omat = $gene.'_tissue.mat';
    open MAT, ">$omat"; print MAT "Gene";
    foreach my $net (@nets) {
        print MAT "\t$net";
    }
    print MAT "\n";

    @other_genes = sort keys %{$gene_val{$gene}};
    foreach my $gene2 (@other_genes) {
        print MAT "$gene2";
        foreach my $val (@{$gene_val{$gene}{$gene2}}) {
            print MAT "\t$val";
        }
        print MAT "\n";
    }

    close MAT;
}

$time = runtime(); print "\n$time: DONE\n\n";

__END__

=head1

Calculate network neighborhood across a number of networks

=head1 USAGE

aggregate_network-neighborhood.pl [--iglist GENES_LIST] [--ibg ALLGENES_LIST]
[--inet NETWORKS_DAB] [--help]

=head1 DESCRIPTION

This script takes in a list of genes & calcualtes its neighborhood in the given
set of networks. Results are aggregated per gene in the form of a matrix with
all network genes along the rows, all the networks along the columns & each cell
(i, j) denoting the connectivity of the gene in question to the row gene i &
network j. The value can additionally be row-z transformed.

=head1 ARGUMENTS

=over 12

=item C<--iglist>

Input list of genes.

=item C<--ibg>

List of all genes (which can be the list of genes in the networks).

=item C<--inet>

Netowrks in DAT/DAB format. Can give more than one network using wild-card
matching.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Jul 12

=cut

