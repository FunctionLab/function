#!/usr/bin/perl

# ================================================
# Name : get_gold-std-genes.pl
# Purpose : Generate a gene-level gold-std for genesets
# Created : 02-01-2014
# Last Modified : Wed 19 Mar 2014 10:20:27 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my ($help, $igmt, $ifiltgs, $iunivg, $irand, $iprior, $inegn, $idir) = @ARGV;
my $irand = 0; my $iming = 5; my $imaxg = 100;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'igmt=s' => \$igmt,
          'idir=s' => \$idir,    
        'iunivg=s' => \$iunivg,
       'ifiltgs=s' => \$ifiltgs,
        'imings=s' => \$imings,
        'imaxgs=s' => \$imaxgs,
         'irand=i' => \$irand,
        'iprior=s' => \$iprior,
         'inegn=i' => \$inegn   ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);
pod2usage( -exitstatus => 2, -verbose => 2 ) unless ($igmt and $idir);

$time = runtime(); print "\n$time: Parsing genes and genesets ...";

# Parse selected genes
my %sel_genes = ();
if($iunivg) {
    %sel_genes = get_selec($iunivg); }

# Parse GMT
my @ref = parse_gmt($igmt, \%sel_genes);
my %gs_posg = %{$ref[0]};
my %gs_desc = %{$ref[1]};
my %gs_size = %{$ref[2]};

$time = runtime(); print "\n$time: Obtaining +ves and –ves ...";

# Assign true gene-labels and record all genes
my %gs_gene_truel = my %all_genes = ();
foreach my $gs (keys %gs_posg) {
    foreach my $g (keys %{$gs_posg{$gs}}) {
        $gs_gene_truel{$gs}{$g} = 1;
        $all_genes{$g}++; } }

# Record slim-term to gene mapping
my %gs_slim = my %slim_genes = ();
if($islim) {
    %gs_slim = %{$ref[4]};
    foreach my $gs (keys %gs_slim) {
        foreach my $s (keys %{$gs_slim{$gs}}) {
            foreach my $g (keys %{$gs_posg{$gs}}) {
                $slim_genes{$s}{$g}++; } } } }

# Parse selected genesets
my %sel_gs = ();
if($igsf) {
    %sel_gs = get_selec($igsf); }

# Get positives and negatives for each geneset
my %gs_negg = my %gs_negg_cvint = ();
my %genes_to_avoid;
foreach my $gs (keys %gs_posg) {
    if($igsf) { unless(exists $sel_gs{$gs}) { next; } }
    %genes_to_avoid = ();

    if($islim) {
        foreach my $s (keys %{$gs_slim{$gs}}) {
            foreach my $g (keys %{$slim_genes{$s}}) {
                $genes_to_avoid{$g}++; } } }

    foreach my $g (keys %all_genes) {
        if(exists $genes_to_avoid{$g}) { next; }
        if(exists $gs_posg{$gs}{$g}) { next; }
        $gs_gene_truel{$gs}{$g} = -1;
        $gs_negg{$gs}{$g}++; }

    $gs_negg_cvint{$gs} = get_cvint($gs_negg{$gs}, $icvk); }



__END__

=head1

Generate gene-level gold-std for genesets.

=head1 USAGE
./get_gold-std-genes.pl [--igmt GENESETS_GMT] [--idir OUTPUT_DIRECTORY]
[--iunivg UNIV_GENES] [--ifiltgs FILTER_GENESETS] [--iming MIN_NUM_GENES]
[--imaxg MIN_NUM_GENES] [--irand] [--iprior PRIOR] [--inegn NUM_NEG_GENES]
[--help]

=head1 DESCRIPTION

This script takes in a collection of genesets in GMT format and outputs a
directory of files, each containing genes labelled as +1, -1 and 0 based on
whether the gene is positive, negative and unknown w.r.t. to each geneset.

=head1 ARGUMENTS

=over 12

=item C<--igmt>

Input collection of genesets in GMT format.

=item C<--idir>

Directory of labelled genes per geneset.

=item C<--iunivg>

Universe of genes in the genome.

=item C<--ifiltgs>

(Optional) List of genesets to filter.

=item C<--imings>

(Optional) Genesets with less than --imings genes will be ignored. Default 5.

=item C<--imaxgs>

(Optional) Genesets with more than --imaxgs genes will be ignored. Default 100.

=item C<--irand>

(Optional) Parameter to set negatives: '0' for the choosing negatives based on
the ontology or slim-terms, '1' for choosing random negatives. Default '0'.

=item C<--iprior>

(Optional) Prior probability of positives to negatives. Default 0.05.

=item C<--inegn>

(Optional) Fixed number of negatives, alternative to choosing a prior.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2014 Jan 02

=cut

