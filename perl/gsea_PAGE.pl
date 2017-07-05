#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw( runtime runinterval figuretimes );
use List::Util 'shuffle';
use Statistics::Distributions;
use Scalar::Util qw(looks_like_number);

sub get_idxarray;
sub parse_gmt;
sub get_gsscore;
sub get_esfdr;
sub clean_fn;

my ($help, $imat, $igmt, $ibg, $ibgco, $omat, $ozmat, $oqmat, $ogenes, $iabs, $itail, $iqval, $iqmat);
my $iming = 10; my $imaxg = 200; my $icol = 'all';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 3);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'icol=s' => \$icol,
          'igmt=s' => \$igmt,
           'ibg=s' => \$ibg,
           'ibgco' => \$ibgco,
         'iming=i' => \$iming,
         'imaxg=i' => \$imaxg,
            'iabs' => \$iabs,
         'itail=s' => \$itail,
         'iqval=f' => \$iqval,
           'iqmat' => \$iqmat,
          'omat=s' => \$omat,
        'ogenes=s' => \$ogenes) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if($iming < 10) { print "\nMinimum geneset size cannot be < 10.\n\n"; exit; }

# Reading in options and all files
my $time = runtime();
print "\n$time: Reading in options and files ...";

open FH,"$imat"; chomp(my @mat=<FH>); close FH;

my $otag = ''; if($iabs) { $otag = '.abs'; }
if($omat) {
    $ozmat = $omat.$otag.'.zscore.mat'; }
else {
    $ozmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.zscore.mat'; }
open ZMAT,">$ozmat";

if($iqmat) {
    if($omat) {
        $oqmat = $omat.$otag.'.esfdr.mat'; }
    else {
        $oqmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.esfdr.mat'; }
    open QMAT,">$oqmat"; }

my ($idsref, $idxref) = get_idxarray($icol, shift(@mat));
my @all_col = @{$idsref}; my @sel_colidx = @{$idxref};

my $numcol = scalar(@sel_colidx);
print "\n\tInput file: $imat\n";
print "\tNo. rows: $#mat\tNo. cols: $numcol";

# Populating genes scores and background population
$time = runtime();
print "\n$time: Indexing input genes, genesets desc and background lists ...\n";

print ZMAT "gs.id\tgs.desc\tgs.ngenes";
if($iqmat) { print QMAT "gs.id\tgs.desc\tgs.ngenes"; }

my %num_genes = (); my %mean_gene_score = (); my %sumsq_gene_score = ();
foreach my $col (@sel_colidx) {
    $num_genes{$col} = 0;
    $mean_gene_score{$col} = 0;
    $sumsq_gene_score{$col} = 0;

    print ZMAT "\t$all_col[$col]";
    if($iqmat) { print QMAT "\t$all_col[$col]"; } }

if($iqmat) { print QMAT "\n"; }
print ZMAT "\n";


# Indexing background genes
my %gmt_bgg = ();
if($ibg) {
    unless($ibgco) {
        open BG, "$ibg" or die "Can't open $ibg!";
        while (<BG>) {
            if($_ =~ /^#/) { next; }
            chomp($_); $gmt_bgg{$_}++; }
        close BG; } }


my %bg_genes = ();
foreach (@mat) {
    if($_ =~ /^#/) { next; }
    my @p = split '\t', $_;
    $bg_genes{$p[0]}++; }


# Indexing gs descriptions and member genes
$time = runtime(); print "$time: Indexing gs descriptions and genes ...\n";

my @gsref = parse_gmt($igmt, \%bg_genes);
my %gs_genes = %{$gsref[0]};
my %gs_desc = %{$gsref[1]};
my %gs_size = %{$gsref[2]};

if($ibgco) {
    unless($ibg) {
        foreach my $gs (keys %gs_genes) {
            if(($gs_size{$gs} < $iming) or ($gs_size{$gs} > $imaxg)) { next; }
            foreach my $g (@{$gs_genes{$gs}}) {
                $gmt_bgg{$g}++; } } } }


# Mean and SumSq calculation using to Welford's algorithm
my %gene_score = (); my (@p, $old_mean);
foreach (@mat) {
    if($_ =~ /^#/) { next; }
    @p = split '\t', $_;
    if($ibg or $ibgco) {
        unless(exists $gmt_bgg{$p[0]}) {
            delete $bg_genes{$p[0]};
            next; } }

    foreach my $col (@sel_colidx) {
        if($p[$col] =~ /(inf|nan|NaN|NA|Inf)/) { $p[$col] = 'NA'; }
        # unless(looks_like_number($p[$col])) { $p[$col] = 'NA'; }
        if($iabs) { unless($p[$col] eq 'NA') { $p[$col] = abs($p[$col]); } }
        $gene_score{$col}{$p[0]} = $p[$col];

        unless($p[$col] eq 'NA') {
            $num_genes{$col}++;
            
            $old_mean = $mean_gene_score{$col};
            $mean_gene_score{$col} += ($p[$col] - $old_mean) / $num_genes{$col};
            $sumsq_gene_score{$col} +=
                ($p[$col] - $old_mean)*($p[$col] - $mean_gene_score{$col}); } } }

# Calculating sd
my %sd_gene_score = my %skip_colidx = (); my $idx = 0;
foreach my $col (@sel_colidx) {
    $sd_gene_score{$col} = sqrt($sumsq_gene_score{$col}/$num_genes{$col});
    if($sd_gene_score{$col} == 0) {
        $skip_colidx{$idx}++; }
    $idx++; }
# print "\t$col\t$mean_gene_score{$col}\t$sd_gene_score{$col}\n"; }

my @temp_colidx = @sel_colidx; @sel_colidx = ();
for(my $i=0; $i<=$#temp_colidx; $i++) {
    if(exists $skip_colidx{$i}) {
        print "\tskipping col $temp_colidx[$i]: $all_col[$temp_colidx[$i]]. no variance.\n";
        next; }
    push(@sel_colidx, $temp_colidx[$i]); }


# Recording geneset scores, and printing zscores
$time = runtime();
print "$time: Calculating gs scores ...", "\n";

my %gs_score = my %gs_modsize = ();
my %gs_zscore = my %gs_parpval = my %gs_zge3 = ();
my (%gs_idx, $par_sd, $zscore, $pvalue); my $count = 0;
foreach my $gs (keys %gs_genes) {
    if(($gs_size{$gs} < $iming) or ($gs_size{$gs} > $imaxg)) { next; }

    $count++; $gs_idx{$gs} = $count;
    @gsref = get_gsscore($gs_genes{$gs}, $gs_size{$gs},
        \%gene_score, \@sel_colidx);

    print ZMAT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    foreach my $col (@sel_colidx) {
        $gs_score{$col}{$gs} = $gsref[0]->{$col};
        $gs_modsize{$col}{$gs} = $gsref[1]->{$col};

        $par_sd = ($sd_gene_score{$col} / sqrt($gs_modsize{$col}{$gs}));
        $zscore = ($gs_score{$col}{$gs} - $mean_gene_score{$col}) / $par_sd;
		print ZMAT "\t", sprintf("%.3f", $zscore);

        $gs_zscore{$col}{$gs} = $zscore;
        if(abs($zscore) >= 3) { $gs_zge3{$gs} = sprintf("%.3f", $zscore); }

        if($iabs) {
            if($itail eq 'ltail') { $pvalue = sprintf("%.5g",
                    (1 - (Statistics::Distributions::uprob($zscore)))); }
            else { $pvalue = sprintf("%.5g",
                    (Statistics::Distributions::uprob($zscore))); } }
        else {
            if($itail) {
                if($itail eq 'utail') { $pvalue = sprintf("%.5g",
                        (Statistics::Distributions::uprob($zscore))); }
                elsif($itail eq 'ltail') { $pvalue = sprintf("%.5g",
                        (1 - (Statistics::Distributions::uprob($zscore)))); } }
            else { $pvalue = sprintf("%.5g",
                    (Statistics::Distributions::uprob(abs($zscore)))*2); } }

        $gs_parpval{$col}{$gs} = $pvalue;
    }
    print ZMAT "\n";
}
close ZMAT;

if($ogenes) {
    my %gene_sym = my %gene_name = ();
    my $igenesym = '/home/arjunk/data/mappings/human_genes_id-symbol-description.txt';
    open ID, "$igenesym" or die "Can't open $igenesym!";
    while (<ID>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $gene_sym{$p[0]} = $p[1];
        $gene_name{$p[0]} = $p[2]; }
    close ID;

    open OG, ">$ogenes";
    print OG "#gs.id\tgs.desc\tgs.ngenes\tgs.zscore\tgene\tsymbol\tscore\tname\n";
    foreach my $gs (sort {$gs_zge3{$b} <=> $gs_zge3{$a}} keys %gs_zge3) {
        my $col = $sel_colidx[0];
        
        my @genea = ();
        if($gs_zge3{$gs} > 0) {
            @genea = sort {$gene_score{$col}{$b} <=> $gene_score{$col}{$a}} @{$gs_genes{$gs}}; }
        else {
            @genea = sort {$gene_score{$col}{$a} <=> $gene_score{$col}{$b}} @{$gs_genes{$gs}}; }

        foreach my $g (@genea) {
            my $y = my $n = '--';
            if(exists $gene_sym{$g}) { $y = $gene_sym{$g}; $n = $gene_name{$g}; }
            print OG "$gs\t$gs_desc{$gs}\t$gs_size{$gs}\t$gs_zge3{$gs}\t$g\t$y\t$gene_score{$col}{$g}\t$n\n"; } }
    close OG; }

my %gs_paresfdr = ();
if($iqval) {
    # Multiple Hypotheses Testing correction
    $time = runtime(); print "$time: Multiple hypothesis correction ...\n";

    my %gs_sigcount = ();
    foreach my $col (@sel_colidx) {
        $gs_paresfdr{$col} = get_esfdr($gs_parpval{$col});

        foreach my $gs (keys %{$gs_paresfdr{$col}}) {
            if($gs_paresfdr{$col}{$gs} >= (-1)*log($iqval)/log(10)) {
                $gs_sigcount{$gs}++; } } }

    # Printing results
    $time = runtime(); print "$time: Printing trimmed Z-score matrix ...\n";

    (my $otzmat = $ozmat) =~ s/\.zscore/\.q$iqval\.zscore/g;
    open TZMAT, ">$otzmat";

    open ZMAT, "$ozmat" or die "Can't open $ozmat!";
    while (<ZMAT>) {
        if($_ =~ /^gs\.id/) { print TZMAT "$_"; next; }
        chomp($_); @p = split '\t', $_;
        unless(exists $gs_sigcount{$p[0]}) { next; }
        if($gs_sigcount{$p[0]} == 0) { next; }

        print TZMAT "$_\n"; }
    close ZMAT; `rm -f $ozmat`;

    close TZMAT; }

if($iqmat) {
    $time = runtime(); print "$time: Printing ES-FDR matrix ...\n";

    my $sign;
    foreach my $gs (sort {$gs_idx{$a} <=> $gs_idx{$b}} keys %gs_idx) {
        print QMAT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
        foreach my $col (@sel_colidx) {
            $sign = 1;
            unless($gs_score{$col}{$gs} == 0) {
                $sign = ($gs_score{$col}{$gs} / abs($gs_score{$col}{$gs})); }
            print QMAT "\t", $sign*$gs_paresfdr{$col}{$gs}; }
        print QMAT "\n"; }
    close QMAT; }

$time = runtime(); print "$time: DONE ...\n\n";

# Subroutines
# ============
# Creates array of indices from a string representation
sub get_idxarray {
    my $string = shift;
    my $header = shift;
    
    my @all_ids = split '\t', $header;
    my @sel_idx = (); my (@p, @q);
    
    if($string eq 'all') {
        for(my $j=1;$j<=$#all_ids;$j++) {
            push(@sel_idx, $j); } }
    else {
        @p = split '_', $string;
        foreach(@p) {
            @q = split ':', $_;
            push(@sel_idx, $q[0]);
            if($#q eq 0) { next; }

            while($q[0]<$q[1]) {
                $q[0]++; push(@sel_idx, $q[0]); } } }

    return (\@all_ids, \@sel_idx); }

# Creates a gs -> genes hash from a .gmt file
sub parse_gmt {
    my $gmt = shift;
    my $bgref = shift;
    
    my (%desc, %genes, %size, @p, $gs);

    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;

        $gs = shift(@p);
        $p[0] =~ s/ \([0-9]*\)$//g;
        $desc{$gs} = $p[0];
        shift(@p);

        foreach my $g (@p) {
            unless(exists $bgref->{$g}) { next; }
            push(@{$genes{$gs}}, $g);
            $size{$gs}++; } }
    close GMT;

    return (\%genes, \%desc, \%size); }

# Calclates mean gs score in all chosen columns
sub get_gsscore {
    my $generef = shift;
    my $size = shift;
    my $scoreref = shift;
    my $colref = shift;
    
    my %mean_score=(); my %mod_size=();
    my $totscore = 0;

    foreach my $col (@{$colref}) {
        $totscore = 0;
        $mod_size{$col} = $size;

        foreach my $gene (@{$generef}) {
            if($scoreref->{$col}->{$gene} eq 'NA') {
                $mod_size{$col}--; }
            else { $totscore += $scoreref->{$col}->{$gene}; } }

        if($mod_size{$col} < $iming) {
            $mean_score{$col} = $mean_gene_score{$col}; }
        else {
            $mean_score{$col} = ($totscore/$mod_size{$col}); } }
    
    return (\%mean_score, \%mod_size); }

# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $pvalues_ref = shift;
    
    my $min_pval = 1; my $num_tests = 0;
    foreach my $gs (keys %{$pvalues_ref}) {
        $num_tests++;
        if(($pvalues_ref->{$gs} != 0) 
                and ($pvalues_ref->{$gs} < $min_pval)) {
            $min_pval = $pvalues_ref->{$gs}; } }
    # print "\tcorrecting for $num_tests\n";

    my %esfdr = (); my $rank = $num_tests;
    foreach my $gs (sort {$pvalues_ref->{$b} <=> $pvalues_ref->{$a}}
        keys %{$pvalues_ref}) {
        if($pvalues_ref->{$gs} == 0) {
            $pvalues_ref->{$gs} = 0.1*sprintf("%.6g", $min_pval); }

        $esfdr{$gs} = sprintf("%.3f",
            (-1)*log(($pvalues_ref->{$gs} * $num_tests) / $rank)/log(10));

        if($esfdr{$gs} < 0) { $esfdr{$gs} = 0; }

        $rank--; }

    return \%esfdr; }

# Cleans filename
sub clean_fn {
    my $str = shift;
    $str =~ s/\.\.\///g;
    $str =~ s/^.*\///g;
    $str =~ s/\.\w*$//g;
    return $str; }


__END__

=head1

Performs Parametric Analysis of Geneset Enrichment (PAGE) based on aggregate of
member-gene scores.

=head1 USAGE

./gsea_PAGE.pl [--imat INPUT_MAT] [--icol COLUMNS] [--igmt GENESETS_GMT] [--help]

=head1 DESCRIPTION

This script takes a matrix containing 'scores' for all profiled genes (rows)
under any number of different factors (columns), and assesses how unlikely are
the aggregate scores for genes in each geneset in the provided collection.

=head1 ARGUMENTS

=over 12

=item C<--imat>

Input matrix. Contains one header row and the first column contains gene ids.
Can contain NAs.

=item C<--igmt>

Collection of genesets in GMT format.

=item C<--icol>

(Optional) String to specify columns of interest in the input matrix using ':'
to specific a continuous series & '_' to specify a break. Index begins at 0.
E.g. "3_5:8_11_13" would mean columns 3, 5, 6, 7, 8, 11 & 13. Default is to use
all columns.

=item C<--iming>

(Optional) Min. no. of genes in a geneset. Genesets with less than <iming> genes
will be excluded. Cannot be < 10. Default 10.

=item C<--imaxg>

(Optional) Max. no. of genes in a geneset. Genesets with greater than <imaxg>
genes will be excluded. Default 200.

=item C<--iabs>

(Optional) Take absolute values of the gene scores & carry-out a one-sided test.

=item C<--itail>

(Optional) '2tail', 'utail' or 'ltail' for calculating two-tail, upper or lower
tail p-value. Default '2tail'.

=item C<--iqval>

(Optional) Q-value cut-off. If provided, a trimmed z-score matrix (removing
genesets that are not significant for any of the columns at provided q-value
cutoff) are printed.

=item C<--omat>

(Optional) Output file tag,  used to create the output filename.
<tag>.zscore.mat will contain the z-scores of genesets along the rows for each
selected input column.
<tag>.q<--iqval>.zscore.mat will contain z-scores only for genesets that have a
q-values < $iqval in at least one column.

=item C<--iqmat>

(Optional) An additional matrix <tag>.esfdr.mat that contains the -log(q-value) of
the genesets for each input column will also be generated.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2011 Feb 15

=cut

