#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Graph::Undirected;
use Time::SoFar qw(runtime);
use PDF;
use Statistics::Distributions;

sub parse_net;
sub parse_gmt;
sub filter_net;
sub get_esfdr;

my ($help, $inet, $igmt1, $igmt2, $iexgs, $otab);
my $imeth = 'bpln'; my $imin = 5; my $imax = 300; my $iqval = 0.01;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'inet=s' => \$inet,
         'igmt1=s' => \$igmt1,
         'igmt2=s' => \$igmt2,
          'imin=i' => \$imin,
          'imax=i' => \$imax,
         'iexgs=s' => \$iexgs,
         'imeth=s' => \$imeth,
         'iqval=i' => \$iqval,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $time, @ref);

# parsing network DAT file
my %net_genes = parse_net($inet);
my $ntotg = scalar keys %net_genes;

# indexing node degree & binning nodes by degree
my %gene_deg = ();
foreach my $g (keys %net_genes) {
    $gene_deg{$g} = scalar keys %{$net_genes{$g}}; }

my $nbin = 20; my $binsize = int($ntotg/$nbin + 0.5);
my $bin = my $idx = 1; my %deg_genes = ();
foreach my $g (sort {$gene_deg{$a} <=> $gene_deg{$b}} keys %gene_deg) {
    if(($idx <= $binsize) or ($bin == $nbin)) { $deg_genes{$bin}{$g}++; }
    else { $binsize += $binsize; $bin++; $deg_genes{$bin}{$g}++; }
    $idx++; }

# parsing genesets GMT files
@ref = parse_gmt($igmt1, \%net_genes);
my %gs1_genes = %{$ref[0]};
my %gs1_desc = %{$ref[1]};
my %gs1_size = %{$ref[2]};

unless($igmt2) { $igmt2 = $igmt1; }
@ref = parse_gmt($igmt2, \%net_genes);
my %gs2_genes = %{$ref[0]};
my %gs2_desc = %{$ref[1]};
my %gs2_size = %{$ref[2]};

%net_genes = filter_net(\%net_genes, \%gs1_genes);

my %exgs_pairs = ();
if($iexgs) {
# indexing geneset pairs to be excluded
    open EX, "$iexgs" or die "Can't open $iexgs!";
    while (<EX>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $exgs_pairs{join '__', sort @p}++;
    }
    close EX;
}

# calculating geneset associations
my $nrand = 100;
my (%ngbh_gs1, %ngbh_gs2, $nngbh_gs1, $nngbh_gs2, $nasce, $tag);
my ($djgs1_size, $djgs2_size, $pval, $zij, $zji);
my %gspair_nasce = (); my %gspair_zscore = ();

GS1: foreach my $gs1 (keys %gs1_genes) {
    GS2: foreach my $gs2 (keys %gs2_genes) {
        if($gs1 eq $gs2) { next; }
        $tag = join '__', sort ($gs1, $gs2);
        if(exists $exgs_pairs{$tag}) { next; }

        if($imeth eq 'bpln') {
            %ngbh_gs1 = (); %ngbh_gs2 = ();
            $djgs1_size = $djgs2_size = 0;

            foreach my $g (keys %{$gs1_genes{$gs1}}) {
                if(exists $gs2_genes{$gs2}{$g}) { next; }
                $djgs1_size++; # gs1 genes not annotated to gs2

                foreach my $ng (keys %{$net_genes{$g}}) {
                    if(exists $gs1_genes{$gs1}{$ng}) { next; }
                    $ngbh_gs1{$ng}++; # net-neighbor of gs1 genes not annotated to gs1
                }
            }

            foreach my $g (keys %{$gs2_genes{$gs2}}) {
                if(exists $gs1_genes{$gs1}{$g}) { next; }
                $djgs2_size++; # gs2 genes not annotated to gs1

                foreach my $ng (keys %{$net_genes{$g}}) {
                    if(exists $gs2_genes{$gs2}{$ng}) { next; }
                    $ngbh_gs2{$ng}++; # net-neighbor of gs2 genes not annotated to gs2
                }
            }

            $nngbh_gs1 = scalar keys %ngbh_gs1;
            $nngbh_gs2 = scalar keys %ngbh_gs2;

            $nasce = 0;
            foreach (keys %ngbh_gs1) {
                if( exists $gs2_genes{$gs2}{$_} ) {
                    $nasce++; } }

            if($nasce < 3) { next GS2; }
            
            $pval = sprintf("%.6g", hypergeometric_tail($ntotg, $nngbh_gs1, $djgs2_size, $nasce));
            if(($pval <= 0) or ($pval >= 1)) { print "ij: $gs1 $gs2: $ntotg, $nngbh_gs1, $djgs2_size, $nasce\n"; exit; }
            $zij = Statistics::Distributions::udistr( hypergeometric_tail($ntotg, $nngbh_gs1, $djgs2_size, $nasce) );

            $nasce = 0;
            foreach (keys %ngbh_gs2) {
                if( exists $gs1_genes{$gs1}{$_} ) {
                    $nasce++; } }

            if($nasce < 3) { next GS2; }

            $pval = sprintf("%.6g", hypergeometric_tail($ntotg, $nngbh_gs2, $djgs1_size, $nasce));
            if(($pval == 0) or ($pval >= 1)) { print "ji: $gs1 $gs2: $ntotg, $nngbh_gs2, $djgs1_size, $nasce\n"; exit; }
            $zji = Statistics::Distributions::udistr( hypergeometric_tail($ntotg, $nngbh_gs2, $djgs1_size, $nasce) );

            $gspair_zscore{$tag} = ($zij + $zji)/sqrt(2);

            $nasce = 0;
            foreach my $g1 (keys %{$gs1_genes{$gs1}}) {
                if(exists $gs2_genes{$gs2}{$g1}) { next; }
                foreach my $g2 (keys %{$gs2_genes{$gs2}}) {
                    if(exists $gs1_genes{$gs1}{$g2}) { next; }
                    if(exists $net_genes{$g1}{$g2}) {
                        $nasce++; } } }
            $gspair_nasce{$tag} = $nasce;
        }

        elsif($imeth eq 'pcna') {
            $nasce = 0;
            foreach my $g1 (keys %{$gs1_genes{$gs1}}) {
                if(exists $gs2_genes{$gs2}{$g1}) { next; }
                foreach my $g2 (keys %{$gs2_genes{$gs2}}) {
                    if(exists $gs1_genes{$gs1}{$g2}) { next; }
                    if(exists $net_genes{$g1}{$g2}) { $nasce++; }
                }
            }

            $gspair_zscore{$tag} = $nasce;
        }
    }
}

open OUT, ">$otab";
print OUT "#GS1.ID\tGS1.Size\tGS1.Desc\tGS2.ID\tGS2.Size\tGS2.Desc\t";
print OUT "No.Edges\tZ.Score\tQ.value\n";

my %gspair_esq = get_esfdr(\%gspair_zscore);
my ($gs1, $gs2); $iqval = (-1)*log($iqval)/log(10);
foreach ( sort {$gspair_esq{$b} <=> $gspair_esq{$a}} keys %gspair_esq ) {
    if($gspair_esq{$_} < $iqval) { last; }
    ($gs1, $gs2) = split '__', $_;

    print OUT "$gs1\t$gs1_size{$gs1}\t$gs1_desc{$gs1}\t";
    print OUT "$gs2\t$gs2_size{$gs2}\t$gs2_desc{$gs2}\t";
    print OUT "$gspair_nasce{$_}\t$gspair_zscore{$_}\t$gspair_esq{$_}\n";
}

close OUT;


# Subroutines
# Get nodes and edges from network
sub parse_net {
    my $dat = shift;

    # my $g = Graph::Undirected->new;
    my %genes = (); my @q;

    open DAT, "$dat" or die "Can't open $dat!";
    while(<DAT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @q = split '\t', $_;
        if($q[0] eq $q[1]) { next; } # skip self-edges
        if($q[2] > 0) {
            # $g->add_edge($q[0], $q[1]);
            $genes{$q[0]}{$q[1]}++; $genes{$q[1]}{$q[0]}++; } }
    close DAT;

    # return ($g, \%genes);
    return %genes; }

# Assigns genes, description & size to gs
sub parse_gmt {
    my $gmt = shift;
    my $bgref = shift;

    my %genes = (); my %desc = (); my %size = ();
    my (@q, $gs);

    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        @q = split '\t', $_;
        $gs = shift(@q);
        ($desc{$gs} = shift(@q)) =~ s/ \([0-9]*\)//g;

        foreach my $g (@q) {
            unless(exists $bgref->{$g}) { next; }
            $genes{$gs}{$g}++; }

        if(scalar keys %{$genes{$gs}} == 0) {
            delete $desc{$gs}; next; }

        $size{$gs} = scalar keys %{$genes{$gs}};
        if(($size{$gs} < $imin) or ($size{$gs} > $imax)) { delete $genes{$gs}; } }
    close GMT;

    return (\%genes, \%desc, \%size); }

# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $zref = shift;
    my $ntests = scalar keys %{$zref};
    my %esfdr = ();

    my $rank = $ntests; my $pval;
    foreach (sort {$zref->{$a} <=> $zref->{$b}} keys %{$zref}) {
        $pval = Statistics::Distributions::uprob($zref->{$_});
        $esfdr{$_} = (-1)*log(($pval * $ntests) / $rank)/log(10);
        if($esfdr{$_} < 0) { $esfdr{$_} = 0; }

        $rank--; }

    return %esfdr; }


__END__

=head1

Calculate associations between genesets in the context of an unweighted network.

=head1 USAGE

calculate_geneset-network-associations.pl [--inet NETWORK_DAT] [--igmt1
GENESETS1_GMT] [--igmt2 GENESETS2_GMT] [--imin MIN_NUM_GENES] [--imax
MAX_NUM_GENES] [--iexgs GENESET_PAIRS_TO_EXCL] [--imeth METHOD] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in two collections of genesets & an unweighted network file to
calculate iassociations between all pairs of genesets across the collections.

=head1 ARGUMENTS

=over 12

=item C<--inet>

Input network in DAT format. Currently only deals with unweighted networks.
Therefore, any edge with non-zero weight is treated as a valid edge (1) & all
othe edges, including gene pairs not explicitly included in the file will be
treated as absent edges (0).

=item C<--igmt1>

Genesets file1 in GMT format.

=item C<--igmt2>

Genesets file2 in GMT format.

=item C<--imin>

(Optional) Gensesets with #genes < imin will be excluded from the analsysis.
Default: 5.

=item C<--imax>

(Optional) Gensesets with #genes > imax will be excluded from the analsysis.
Default: 300.

=item C<--iexgs>

(Optional) List of geneset pairs to avoid. This could be the list of geneset
pairs that indeed have a significant number of overlapping genes that the
association can be explained in terms of gene-overlap instead of 'cross-talk'.

=item C<--imeth>

Method for calculating geneset associations.

[bpln] Parametric method from Dotan-Cohen (2009) PLoS ONE 4: e5313.

[pcna] Non-parametric method modified from Li (2008) Bioinformatics 24: 1442.

=item C<--iqval>

(Optional) Q-value (FDR) cut-off to declare a geneset pair significant. Default
is 0.01.

=item C<--otab>

Output table reporting all statistics and similarity scores.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Jul 03

=cut

