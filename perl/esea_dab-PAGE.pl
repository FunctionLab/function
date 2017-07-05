#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Time::SoFar qw( runtime );
use Statistics::Distributions;

my ($help, $idab, $igmt, $iskip, $icost, $icover, $otab);
my $itail = '2tail'; my $imiss = 'NA'; my $imode = 'both';
my $icdir = './contexts'; my $iodir = '.';
my $iming = 5; my $imaxg = 200; my $imine = 10; my $imaxe = 500;

pod2usage( -exitstatus => 2, -verbose => 2 ) if ( @ARGV < 2 );
GetOptions( 'help' => \$help,
          'idab=s' => \$idab,
          'igmt=s' => \$igmt,
         'iming=s' => \$iming,
         'imaxg=s' => \$imaxg,
         'imine=s' => \$imine,
         'imaxe=s' => \$imaxe,
         'icdir=s' => \$icdir,
         'iskip=s' => \$iskip,
         'icost=s' => \$icost,
        'icover=s' => \$icover,
         'itail=s' => \$itail,
         'imiss=f' => \$imiss,
         'imode=s' => \$imode,
         'iodir=s' => \$iodir,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $time); $icdir =~ s/\/$//g; $iodir =~ s/\/$//g;


$time = runtime();
print "\n$time: Indexing gs genes and desc ...\n";

my %gs_genes = my %gs_desc = ();

open GMT, "$igmt" or die "Can't open $igmt!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    my $gs = shift(@p);
    ($gs_desc{$gs} = shift @p) =~ s/ \([0-9]*\)$//g;
    $gs_desc{$gs} =~ s/ /_/g; $gs_desc{$gs} = lc($gs_desc{$gs});

    foreach my $g (@p) {
        $gs_genes{$gs}{$g}++; } }
close GMT;


$time = runtime();
print "\n$time: Indexing input edges, scores and background lists ...\n";

(my $igenes = $idab) =~ s/\.dab/\.genes/g; $igenes =~ s/^.*\///g;
$igenes = $iodir. '/' . $igenes;
`Dat2Dab -i $idab -P > $igenes`;

open GEN, "$igenes" or die "Can't open $igenes!";
chomp(my @g = <GEN>); close GEN;

my $ntote = $g[0];
my ($mean_escore, $sd_escore) = split '\t', $g[1];
shift @g; shift @g;

my %net_genes = ();
foreach (@g) {
    @p = split '\t', $_;
    $net_genes{$p[0]} = $p[1]; }

print "\tTotal no. genes: ", scalar keys %net_genes, "\n";
print "\tTotal no. edges: $ntote\n";
print "\tMean escore: $mean_escore\n\tSD escores: $sd_escore\n";


$time = runtime();
print "\n$time: Printing contexts and indexing co-ann edges ...\n";

unless(-d $icdir) { `mkdir $icdir/`; }
unless(-d $iodir) { `mkdir $iodir/`; }

my %gs_gsize = my %gs_esize = my %gs_escore = ();
my @gs_array = my @final_gs_array = ();

foreach my $gs (keys %gs_genes) {
    $gs_gsize{$gs} = 0;
    foreach my $g (keys %{$gs_genes{$gs}}) {
        if(exists $net_genes{$g}) { $gs_gsize{$gs}++; }
        else { delete $gs_genes{$gs}{$g}; } }

    if(($gs_gsize{$gs} < $iming) or ($gs_gsize{$gs} > $imaxg)) {
        delete $gs_genes{$gs}; next; }
    push(@gs_array, $gs);
    
    $time = runtime();
    print "\t$time: $gs_desc{$gs}\n";

    my @tgenes = keys %{$gs_genes{$gs}};

    my $tdesc = $gs_desc{$gs};
    my $tctxt = $icdir. '/' .$tdesc;

    unless(-e $tctxt) {
        open CXT, ">$tctxt";
        foreach my $g (@tgenes) {
            print CXT "$g\n"; }
        close CXT; }

    if(($imode eq 'intra') || ($imode eq 'both')) {
        my $tdat = $iodir . '/' . $tdesc . '.dat';
        `Dat2Dab -i $idab -g $icdir/$tdesc -o $tdat`;

        $gs_escore{$gs} = $gs_esize{$gs} = 0;
        open DAT, "$tdat" or die "Can't open $tdat!";
        while (<DAT>) {
            if($_ =~ /^#/) { next; }
            chomp($_); @p = split '\t', $_;

            $gs_esize{$gs}++;
            $gs_escore{$gs} += $p[2]; }
        close DAT;

        $gs_escore{$gs} = ($gs_escore{$gs} / $gs_esize{$gs});

        if(($gs_esize{$gs} < $imine) or ($gs_esize{$gs} > $imaxe)) { next; }
        push(@final_gs_array, $gs); } }


if(($imode eq 'inter') || ($imode eq 'both')) {
    $time = runtime();
    print "\n$time: Indexing cross-ann edges ...\n";

    my %gs_toskip = ();
    if($iskip) {
        open SK, "$iskip" or die "Can't open $iskip!";
        while (<SK>) {
            if($_ =~ /^#/) { next; }
            chomp($_); @p = split '\t', $_;

            my $gs = $p[0] . '|' . $p[3];
            $gs_toskip{$gs}++; }
        close SK; }

    for(my $i=0; $i<$#gs_array; $i++) {
        my $gs1 = $gs_array[$i];
        my $tdesc1 = $gs_desc{$gs1};

        for(my $j=($i+1); $j<=$#gs_array; $j++) {
            my $gs2 = $gs_array[$j];
            my $tdesc2 = $gs_desc{$gs2};

            my $gs = $gs1 . '|' . $gs2;
            if($iskip) { if(exists $gs_toskip{$gs}) { next; } }

            my $tdat = $iodir . '/' . $gs_desc{$gs1} . '--' . $gs_desc{$gs2} . '.dat';
            `Dat2Dab -i $idab -t $icdir/$tdesc1 -T $icdir/$tdesc2 > $tdat`;

            my %genes = ();
            $gs_escore{$gs} = $gs_esize{$gs} = 0;
            open DAT, "$tdat" or die "Can't open $tdat!";
            while (<DAT>) {
                if($_ =~ /^#/) { next; }
                chomp($_); @p = split '\t', $_;
                $genes{$p[0]}++; $genes{$p[1]}++;

                $gs_esize{$gs}++;
                $gs_escore{$gs} += $p[2]; }
            close DAT;

            if(($gs_esize{$gs} < $imine) or ($gs_esize{$gs} > $imaxe)) { next; }
            push(@final_gs_array, $gs);

            $gs_desc{$gs} = $tdesc1. '|' . $tdesc2;
            $gs_gsize{$gs} = scalar keys %genes;
            $gs_escore{$gs} = ($gs_escore{$gs} / $gs_esize{$gs}); } } }


$time = runtime(); print "\n$time: Calculating gs scores ...", "\n";

my %gs_mscore = (); my %gs_sscore = ();
my %gs_zscore = (); my %gs_pvalue = ();

foreach my $gs (@final_gs_array) {
	my $sd = ($sd_escore / sqrt($gs_esize{$gs}));

    my $zscore = (($gs_escore{$gs} - $mean_escore) / $sd);
    my $pvalue = get_pvalue($zscore, $itail);

	$gs_mscore{$gs} = sprintf("%.6f", $gs_escore{$gs});
	$gs_sscore{$gs} = sprintf("%.6f", ($gs_escore{$gs}*$gs_esize{$gs}));
	$gs_zscore{$gs} = sprintf("%.6f", $zscore);
	$gs_pvalue{$gs} = sprintf("%.6g", $pvalue); }


$time = runtime(); print "\n$time: Calculating corrected qvalues and printing results ...\n";

my %gs_esfdr = get_esfdr(\%gs_pvalue);

open TAB, ">$iodir/$otab";
print TAB "#GS\tDesc\tNo.Genes\tNo.Edges\tSum\tMean\tZscore\tesFDR\n";
foreach my $gs (sort {$gs_zscore{$b} <=> $gs_zscore{$a}} keys %gs_esfdr) {
	print TAB "$gs\t$gs_desc{$gs}\t$gs_gsize{$gs}\t$gs_esize{$gs}";
	print TAB "\t$gs_sscore{$gs}\t$gs_mscore{$gs}\t$gs_zscore{$gs}\t$gs_esfdr{$gs}\n"; }
close TAB;


$time = runtime();
print "\n$time: DONE\n\n";


# Subroutines
sub get_gsscore {
    my $eref = shift;
    my $size = shift;
    my $esco = shift;

    my $tots = 0;
    foreach my $e (keys %{$eref}) {
        $tots += $esco->{$e}; }

    my $mean = ( $tots / $size );
    return $mean; }

sub get_pvalue {
    my $zsco = shift;
    my $tail = shift;

    my $pval = 1;
	if($tail eq '2tail') {
        $pval = sprintf("%.6g", (Statistics::Distributions::uprob (abs($zsco)))*2); }
	elsif($tail eq 'utail') {
        $pval = sprintf("%.6g", (Statistics::Distributions::uprob ($zsco))); }
	elsif($tail eq 'ltail') {
        $pval = sprintf("%.6g", (1 - Statistics::Distributions::uprob ($zsco))); }

    return $pval; }

# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $pref = shift;
    my $min_pval = 1;
    my $ntest = 0;
    
    foreach my $gs (keys %{$pref}) {
        $ntest++;
        if(($pref->{$gs} != 0) 
                and ($pref->{$gs} < $min_pval)) {
            $min_pval = $pref->{$gs}; } }

    my %esfdr = (); my $rank = $ntest;
    foreach my $gs (sort {$pref->{$b} <=> $pref->{$a}} keys %{$pref}) {
        if($pref->{$gs} == 0) {
            $pref->{$gs} = 0.1*sprintf("%.6g", $min_pval); }

        $esfdr{$gs} = sprintf("%.6f",
            (-1)*log(($pref->{$gs} * $ntest) / $rank)/log(10));

        if($esfdr{$gs} < 0) { $esfdr{$gs} = 0; }

        $rank--; }

    return %esfdr; }


__END__

=head1

Calculates enrichment of gene sets in edge-space based on an underlying gene
network.

=head1 USAGE

esea_dab-PAGE.pl [--idab NETWORK_DAB] [--igmt GENESETS_GMT] [--icdir CONTEXT_DIR] [--iming MIN_NUM_GENES]
[--imaxg MAX_NUM_GENES] [--imine MIN_NUM_EDGES] [--imaxe MAX_NUM_EDGES] [--itail
TAIL_OF_TEST] [--imiss MISSING_VALUES] [--imode ENRICHMENT_MODE] [--iskip
GENESET-PAIRS_TO-SKIP] [--icost GENESET_COSTS] [--icover OUTPUT_GENESET_COVER]
[--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in a gene network and a collection of genesets, and
calculates the edge-weight-based enrichment of intra- and/or inter-geneset
edges.

=head1 ARGUMENTS

=over 12

=item C<--idab>

Weighted network file in DAB format.

=item C<--imiss>

(Optional) Value for missing edges in the network. By default, these edges will be ignored.
If this option is provided, then all missing edges are taken into account and
the value provided will be assigned to these edges as their weights. Not
implemented yet.

=item C<--igmt>

Geneset collection in GMT format.

=item C<--icdir>

(Optional) Directory to print contexts. Default './contexts/'.

=item C<--imode>

(Optional) Edge-types to consider for analysis. 'intra' - consider only edges between genes
that are coannotated to at least one of the selected genesets. 'inter' -
consider only edges between genes that are each annotated to one of the selected
GO terms but not coannotated to any. 'both' - both 'intra' and 'inter' edges.
Default 'both'.

=item C<--iming>

(Optional) Min. no. genes in genesets. Genesets with less genes will be ignored. Default 5.

=item C<--imaxg>

(Optional) Max. no. genes in genesets. Genesets with more genes will be ignored. Default
200.

=item C<--imine>

(Optional) Min. no. edges in genesets. Genesets with less edges will be ignored. Default
10.

=item C<--imaxe>

(Optional) Max. no. edges in genesets. Genesets with less edges will be ignored. Default
500.

=item C<--itail>

(Optional) Tail of the enrichment hypothesis test: '2tail', 'upper' or 'lower'. Default
'2tail'.

=item C<--iskip>

(Optional) Pairs of genesets to skip during an 'inter' analysis mode. These could be
genesets that have a significant gene-overlap which would explain extensive
connectivity between the genesets and not necessarily point to a 'crosstalk'.
This file should be in the following format with each line listing a pair of
genesets: <gs1> <gs1_desc> <gs1_size> <gs2> <gs2_desc> <gs2_size>

=item C<--icost>

(Optional) Tab-delimited file containing costs for each geneset: <GS.ID> <Cost>.
Not implemented yet.

=item C<--icover>

(Optional) Perform set-cover and file to print output. Not implemented yet.

=item C<--iodir>

(Optional) Directory to print output DAT files and tables. Default './'.

=item C<--otab>

Output file containing a table of relevant genesets, their gene and edge
coverage statistics and their enrichment scores.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2010 May 15

=cut

