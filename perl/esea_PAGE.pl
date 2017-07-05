#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Time::SoFar qw( runtime );
use Statistics::Distributions;

my ($help, $idat, $igmt, $iskip, $icost, $icover, $otab);
my $itail = 'twotail'; my $imiss = 'NA'; my $imode = 'both';
my $iming = 5; my $imaxg = 200; my $imine = 10; my $imaxe = 500;

pod2usage( -exitstatus => 2, -verbose => 2 ) if ( @ARGV < 2 );
GetOptions( 'help' => \$help,
          'idat=s' => \$idat,
          'igmt=s' => \$igmt,
         'iming=s' => \$iming,
         'imaxg=s' => \$imaxg,
         'imine=s' => \$imine,
         'imaxe=s' => \$imaxe,
         'iskip=s' => \$iskip,
         'icost=s' => \$icost,
        'icover=s' => \$icover,
         'itail=s' => \$itail,
         'imiss=f' => \$imiss,
         'imode=s' => \$imode,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $e, @tgenes, $gs, $g1, $g2, $time);

$time = runtime();
print "\n$time: Indexing gs genes and desc ...\n";

my %gs_genes = my %gene_gs = my %gs_desc = ();

open GMT, "$igmt" or die "Can't open $igmt!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $gs = shift(@p);
    ($gs_desc{$gs} = shift @p) =~ s/ \([0-9]*\)$//g;

    foreach my $g (@p) {
        $gs_genes{$gs}{$g}++;
        $gene_gs{$g}{$gs}++; } }
close GMT;

$time = runtime();
print "\n$time: Indexing input edges, scores and background lists ...\n";

my %escore = ();
my %net_genes = (); my %net_edges = (); my $ntote = 0;
my $mean_escore = my $old_mean = my $sd_escore = 0;

open DAT, "$idat" or die "Can't open $idat!";
while (<DAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $e = join '__', sort($p[0], $p[1]);
    $ntote++;
    $mean_escore = $old_mean + ($p[2] - $old_mean)/$ntote;
    $sd_escore += ($p[2] - $old_mean)*($p[2] - $mean_escore);
    $old_mean = $mean_escore;

    if((exists $gene_gs{$p[0]}) and (exists $gene_gs{$p[1]})) {
        $net_genes{$p[0]}++; $net_genes{$p[1]}++;
        $net_edges{$e}++; $escore{$e} = $p[2]; } }
close DAT;

my $nmisse = 0;
if($imiss ne 'NA') {
	@tgenes = keys %net_genes;

	for(my $i=0; $i<$#tgenes; $i++) {
		for(my $j=($i+1); $j<=$#tgenes; $j++) {
			$e = join '__', sort($tgenes[$i], $tgenes[$j]);

			unless(exists $net_edges{$e}) {
                $escore{$e} = $imiss; $ntote++;
                $mean_escore = $old_mean + ($imiss - $old_mean)/$ntote;
                $sd_escore += ($imiss - $old_mean)*($imiss - $mean_escore);
                $old_mean = $mean_escore;
				$nmisse++; $net_edges{$e}++; } } } }

$sd_escore = sqrt($sd_escore / $ntote);

print "\tTotal no. edges: $ntote\n";
print "\tMean escore: $mean_escore\n\tSD escores: $sd_escore\n";

$time = runtime();
print "\n$time: Indexing co-ann edges ...\n";

my %gs_gsize = ();
my %gs_edges = (); my %gs_esize = ();
my %ann_edges = (); my %coann_edges = ();
my %gspair_desc = ();

foreach my $gs (keys %gs_genes) {
    $gs_gsize{$gs} = 0;
    foreach my $g (keys %{$gs_genes{$gs}}) {
        if(exists $net_genes{$g}) { $gs_gsize{$gs}++; }
        else { delete $gs_genes{$gs}{$g}; } }

    if(($gs_gsize{$gs} < $iming) or ($gs_gsize{$gs} > $imaxg)) {
        delete $gs_genes{$gs}; next; }

    $gspair_desc{$gs} = $gs_desc{$gs};
    @tgenes = keys %{$gs_genes{$gs}};
	
	for(my $i=0; $i<$#tgenes; $i++) {
		for(my $j=($i+1); $j<=$#tgenes; $j++) {
            $e = join '__', sort($tgenes[$i], $tgenes[$j]);
			unless(exists $net_edges{$e}) { next; }

            $ann_edges{$e}++;
            $coann_edges{$e}{$gs}++;

            if(($imode eq 'intra') or ($imode eq 'both')) {
                $gs_edges{$gs}{$e}++;
                $gs_esize{$gs}++; } } } }

foreach my $g (keys %gene_gs) {
    unless(exists $net_genes{$g}) { delete $gene_gs{$g}; } }

$time = runtime();
print "\n$time: Indexing cross-ann edges ...\n";

my $ncoanne = scalar keys %coann_edges;
my $ncranne = 0;

foreach my $e (keys %net_edges) {
    if(exists $coann_edges{$e}) { next; }
    ($g1, $g2) = split '__', $e;
    unless((exists $gene_gs{$g1}) and (exists $gene_gs{$g2})) { next; }

    $ncranne++;
    unless($ncranne % 100000) {
        $time = runtime(); print "\t$time: $ncranne ...\n"; }

    if($imode eq 'intra') { next; }

    foreach my $gs1 (keys %{$gene_gs{$g1}}) {
        foreach my $gs2 (keys %{$gene_gs{$g2}}) {
            if($gs1 eq $gs2) { next; }
            @p = sort($gs1, $gs2);
            $gs = join '|', @p;
            $gs_edges{$gs}{$e}++;
            $ann_edges{$e}++;

            $gs_esize{$gs}++;
            unless(exists $gspair_desc{$gs}) {
                $gspair_desc{$gs} = $gs_desc{$p[0]}.'|'.$gs_desc{$p[1]}; } } } }

my $nanne = scalar keys %ann_edges;
my $ngs = scalar keys %gs_esize;

$time = runtime();
print "\n$time: Recording intra- and inter-gs counts ...\n";

my $ngs_wenoughe = 0; my @gs_array = ();
my $nintra_gs = my $ninter_gs = 0;

my %dyn_gs_esize = (); my %edge_gs = ();
foreach my $gs (keys %gs_esize) {
    if(($gs_esize{$gs} < $imine) or ($gs_esize{$gs} > $imaxe)) { next; }
    $dyn_gs_esize{$gs} = $gs_esize{$gs};
    $ngs_wenoughe++;

    push(@gs_array, $gs);

    foreach my $e (keys %{$gs_edges{$gs}}) {
        $edge_gs{$e}{$gs}++; }

    if($gs =~ /\|/) {
        $ninter_gs++;
        foreach my $e (keys %{$gs_edges{$gs}}) {
            @p = split '__', $e;
            $gs_genes{$gs}{$p[0]}++; $gs_genes{$gs}{$p[1]}++; }
        $gs_gsize{$gs} = scalar keys %{$gs_genes{$gs}}; }
    else {
        $nintra_gs++; } }

print "\n\tTotal no. ann edges: $nanne\n\tNo. coann edges: $ncoanne\n\tNo. cross-ann edge: $ncranne\n";
print "\n\tTotal no. genesets: $ngs\n\tNo. genesets w/ >= $imine & <= $imaxe edges: $ngs_wenoughe\n";
print "\tNo. intrags: $nintra_gs\n\tNo. intergs: $ninter_gs\n\n";

$time = runtime(); print "$time: Calculating gs scores ...", "\n";

my %gs_mscore = (); my %gs_sscore = ();
my %gs_zscore = (); my %gs_pvalue = ();
my ($mean, $sd, $zscore, $pvalue);

foreach $gs (@gs_array) {
	$sd = ($sd_escore/sqrt($gs_esize{$gs}));

    $mean = get_gsscore($gs_edges{$gs}, $gs_esize{$gs}, \%escore);
    $zscore = (($mean - $mean_escore) / $sd);
    $pvalue = get_pvalue($zscore, $itail);

	$gs_mscore{$gs} = sprintf("%.6f", $mean);
	$gs_sscore{$gs} = sprintf("%.6f", ($mean*$gs_esize{$gs}));
	$gs_zscore{$gs} = sprintf("%.6f", $zscore);
	$gs_pvalue{$gs} = sprintf("%.6g", $pvalue); }

$time = runtime(); print "\n$time: Calculating corrected qvalues ...\n";

my %gs_esfdr = get_esfdr(\%gs_pvalue);

$time = runtime(); print "\n$time: Printing scores ...\n";

open TAB, ">$otab";
print TAB "#GS\tDesc\tNo.Genes\tNo.Edges\tSum\tMean\tZscore\tesFDR\n";
foreach $gs (keys %gs_esfdr) {
	print TAB "$gs\t$gspair_desc{$gs}\t$gs_gsize{$gs}\t$gs_esize{$gs}";
	print TAB "\t$gs_sscore{$gs}\t$gs_mscore{$gs}\t$gs_zscore{$gs}\t$gs_esfdr{$gs}\n"; }
close TAB;

unless($icover) {
	$time = runtime();
	print "$time: DONE\n\n"; exit; }

$time = runtime(); print "\n$time: Performing set-cover ...\n";

my %gs_cost = ();
if($icost) {
    open INF, "$icost" or die "Can't open $icost!";
    while (<INF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        # $gs_cost{$p[0]} = 1 - $p[1]; }
        $gs_cost{$p[0]} = $p[1]; }
    close INF;

    foreach (keys %gs_zscore) {
        if(exists $gs_cost{$_}) { next; }
        unless($_ =~ /\|/) { next; }
        @p = split '\|', $_;
        if((exists $gs_cost{$p[0]}) and ($gs_cost{$p[1]})) {
        $gs_cost{$_} = ($gs_cost{$p[0]} + $gs_cost{$p[1]}) / 2; } } }

open SET, ">$icover";
print SET "#ID\tDesc\tNo.Genes\tNo.Edges\tCost\tSum\tMean\tZscore\tesFDR\t";
print SET "No.New.Edges\tFrac.Cov.NetEdges\tFrac.Cov.GSEdges\tCostBenf.ratio\n";

my %remaining_gs = ();
my %remaining_edges = ();
foreach my $gs (keys %gs_esfdr) {
    if($gs_esfdr{$gs} >= 2) {
        $remaining_gs{$gs}++;
        foreach my $e (keys %{$gs_edges{$gs}}) {
            $remaining_edges{$e}++; } } }

my $num_reme = scalar keys %remaining_edges;
my $num_cove = my $count = 0;
my ($min_costbenf, $max_costbenf, $min_gs, $max_gs, $benf, $costbenf);
my %gs_frace = ();

print "\tTot. no. remaining edges: $num_reme\n\n";

while($num_reme > 0) {
    $count++;
    # $min_costbenf = 1000;
    $max_costbenf = 0;
    foreach my $gs (keys %remaining_gs) {
        # $costbenf = ( log($gs_pvalue{$gs}) / (log(10)*$dyn_gs_esize{$gs}) );
        # $benf = (-1)*log( $dyn_gs_esize{$gs} / $nanne );
        # $costbenf = $gs_cost{$gs} / $benf;
        # $costbenf = $gs_cost{$gs}*log(2) / log($dyn_gs_esize{$gs});
        # $costbenf = $gs_cost{$gs} * log(2) / log($dyn_gs_esize{$gs});
        # $costbenf = ( $gs_cost{$gs} / abs($gs_zscore{$gs}) ) * ( log(2) / log($dyn_gs_esize{$gs}) );
        $costbenf = ($gs_cost{$gs}**2) * abs($gs_zscore{$gs}) * log($dyn_gs_esize{$gs}) / log(2);
        # $costbenf = $gs_cost{$gs} / $dyn_gs_esize{$gs};
        # $costbenf = $gs_cost{$gs}*abs($gs_zscore{$gs}) / $dyn_gs_esize{$gs};
        # if($costbenf < $min_costbenf) {
        if($costbenf > $max_costbenf) {
            # $min_gs = $gs;
            $max_gs = $gs;
            # $min_costbenf = $costbenf; } }
            $max_costbenf = $costbenf; } }

    $num_cove += $dyn_gs_esize{$max_gs};
    $num_reme -= $num_cove;

    if(($dyn_gs_esize{$max_gs} / $gs_esize{$max_gs}) < 0.25) { last; }

    # unless($min_gs) { die "no min_gs!"; }
    # unless($gspair_desc{$min_gs}) { die "no desc min_gs!"; }
    # unless($gs_gsize{$max_gs}) { die "no gsize min_gs!"; }
    # unless($gs_esize{$max_gs}) { die "no esize min_gs!"; }

    print "\t$max_gs\t$gs_esize{$max_gs}\t$dyn_gs_esize{$max_gs}\t$gs_cost{$max_gs}\t$max_costbenf\t$num_cove\t$num_reme\n";

    print SET "$max_gs\t$gspair_desc{$max_gs}\t$gs_gsize{$max_gs}\t$gs_esize{$max_gs}\t$gs_cost{$max_gs}";
    print SET "\t$gs_sscore{$max_gs}\t$gs_mscore{$max_gs}\t$gs_zscore{$max_gs}\t$gs_esfdr{$max_gs}";
    print SET "\t$dyn_gs_esize{$max_gs}", sprintf("\t%.3g\t%.3g", ($dyn_gs_esize{$max_gs}/$gs_esize{$max_gs}), ($num_cove/$nanne));
    print SET sprintf("\t%.3g\n", $max_costbenf);

    foreach my $e (keys %{$gs_edges{$max_gs}}) {
        foreach my $gs (keys %{$edge_gs{$e}}) {
            $dyn_gs_esize{$gs}--; } }

    foreach my $gs (keys %remaining_gs) {
        if($dyn_gs_esize{$gs} < 2) {
            delete $remaining_gs{$gs}; } }

    unless($count % 100) {
        $time = runtime();
        print "\t$time: $max_gs\t$dyn_gs_esize{$max_gs}\t$gs_cost{$max_gs}\t$max_costbenf\t$num_cove\t$num_reme\n"; } }

$time = runtime();
print "\n$time: No. covered / uncovered edges: $num_cove / $num_reme\n\nDONE\n\n";

close SET;


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
	if($tail eq 'twotail') {
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

esea_PAGE.pl [--idat NETWORK_DAT] [--igmt GENESETS_GMT] [--iming MIN_NUM_GENES]
[--imaxg MAX_NUM_GENES] [--imine MIN_NUM_EDGES] [--imaxe MAX_NUM_EDGES] [--itail
TAIL_OF_TEST] [--imiss MISSING_VALUES] [--imode ENRICHMENT_MODE] [--icost
GENESET_COSTS] [--icover OUTPUT_GENESET_COVER] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in a gene network and a collection of genesets, and
calculates the edge-weight-based enrichment of intra- and/or inter-geneset
edges.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Weighted network file in DAT format.

=item C<--imiss>

Value for missing edges in the network. By default, these edges will be ignored.
If this option is provided, then all missing edges are taken into account and
the value provided will be assigned to these edges as their weights.

=item C<--igmt>

Geneset collection in GMT format.

=item C<--imode>

Edge-types to consider for analysis. 'intra' - consider only edges between genes
that are coannotated to at least one of the selected genesets. 'inter' -
consider only edges between genes that are each annotated to one of the selected
GO terms but not coannotated to any. 'both' - both 'intra' and 'inter' edges.
Default 'both'.

=item C<--iming>

Min. no. genes in genesets. Genesets with less genes will be ignored. Default 5.

=item C<--imaxg>

Max. no. genes in genesets. Genesets with more genes will be ignored. Default
200.

=item C<--imine>

Min. no. edges in genesets. Genesets with less edges will be ignored. Default
10.

=item C<--imaxe>

Max. no. edges in genesets. Genesets with less edges will be ignored. Default
500.

=item C<--itail>

Tail of the enrichment hypothesis test: 'twotail', 'upper' or 'lower'. Default
'twotail'.

=item C<--iskip>

Pairs of genesets to skip during an 'inter' analysis mode. These could be
genesets that have a significant gene-overlap which would explain extensive
connectivity between the genesets and not necessarily point to a 'crosstalk'.

=item C<--icost>

Tab-delimited file containing costs for each geneset: <GS.ID>	<Cost>

=item C<--icover>

Option to perform set-cover and file to print output.

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

