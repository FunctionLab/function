#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw( runtime runinterval figuretimes );
use List::Util 'shuffle';
use Statistics::Distributions;

sub get_idxarray;
sub parse_gmt;
sub get_gsscore;
sub get_esfdr;
sub clean_fn;

my ($help, $imat, $igmt, $omat, $ozmat, $oqmat, $otab, $iabs, $itail, $iqval);
my $imin = 10; my $imax = 200; my $icol = 'all';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 3);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'icol=s' => \$icol,
          'igmt=s' => \$igmt,
          'imin=i' => \$imin,
            'iabs' => \$iabs,
         'itail=s' => \$itail,
         'iqval=f' => \$iqval,
          'imax=i' => \$imax,
          'omat=s' => \$omat) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if($imin < 10) { print "\nMinimum geneset size cannot be < 10.\n\n"; exit; }

# Reading in options and all files
my $time = runtime();
print "\n$time: Reading in options and files ...";

open FH,"$imat"; chomp(my @mat=<FH>); close FH;

my $otag = ''; if($iabs) { $otag = '.abs'; }
if($omat) {
    $ozmat = $omat.$otag.'.zscore.mat';
    $oqmat = $omat.$otag.'.esfdr.mat';
    $otab = $omat.$otag.'.zscore-esfdr.txt'; }
else {
    $ozmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.zscore.mat';
    $oqmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.esfdr.mat';
    $otab = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.zscore-esfdr.txt'; }
open ZMAT,">$ozmat"; open QMAT,">$oqmat"; open TAB,">$otab";

# Selecting columns
my %fcol = (); my @p;
open COL, "$icol" or die "Can't open $icol!";
while (<COL>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    $fcol{$p[0]}++; }
close COL;

my $header = shift @mat;
my @all_col = split '\t', $header;

my @sel_colidx = ();
for(my $j=0; $j<=$#all_col; $j++) {
    if(exists $fcol{$all_col[$j]}) {
        push(@sel_colidx, $j); } }
# $sel_colidx{$j} = $all_col[$j]; } }

my $numcol = scalar @sel_colidx;
print "\nInput file: $imat\n";
print "No. rows: ", ($#mat+1), "\tNo. cols: $numcol\n";

# Populating genes scores and background population
$time = runtime();
print "\n$time: Indexing input genes, genesets desc and background lists ...\n";

print ZMAT "GS.ID\tGS.Desc\tGS.Size";
unless($iqval) { print QMAT "GS.ID\tGS.Desc\tGS.Size"; }
print TAB "GS.ID\tGS.Desc\tGS.Size\tZscore\tesFDR\n";

my %num_genes = (); my %mean_gene_score = (); my %sumsq_gene_score = ();
foreach my $col (@sel_colidx) {
    $num_genes{$col} = 0;
    $mean_gene_score{$col} = 0;
    $sumsq_gene_score{$col} = 0;

    print ZMAT "\t$all_col[$col]";
    unless($iqval) { print QMAT "\t$all_col[$col]"; }
}

unless($iqval) { print QMAT "\n"; }
print ZMAT "\n";

# Mean and SumSq calculation using to Welford's algorithm
my %gene_score = (); my %bg_genes = (); my $old_mean;
foreach (@mat) {
    if($_ =~ /^#/) { next; }
    @p = split '\t', $_;
    $bg_genes{$p[0]}++;

    foreach my $col (@sel_colidx) {
        if($iabs) { unless($p[$col] eq 'NA') { $p[$col] = abs($p[$col]); } }
        $gene_score{$col}{$p[0]} = $p[$col];

        unless($p[$col] eq 'NA') {
            $num_genes{$col}++;
            
            $old_mean = $mean_gene_score{$col};
            $mean_gene_score{$col} += ($p[$col] - $old_mean) / $num_genes{$col};
            $sumsq_gene_score{$col} +=
                ($p[$col] - $old_mean)*($p[$col] - $mean_gene_score{$col}); } } }

# Calculating sd
my %sd_gene_score = ();
foreach my $col (@sel_colidx) {
    $sd_gene_score{$col} = sqrt($sumsq_gene_score{$col}/$num_genes{$col}); }

# Indexing gs descriptions and member genes
$time = runtime(); print "$time: Indexing gs descriptions and genes ...\n";

my @gsref = parse_gmt($igmt, \%bg_genes);
my %gs_genes = %{$gsref[0]};
my %gs_desc = %{$gsref[1]};
my %gs_size = %{$gsref[2]};

# Recording geneset scores, and printing zscores
$time = runtime();
print "$time: Calculating gs scores ...", "\n";

my %gs_score = (); my %gs_modsize = ();
my %gs_zscore = (); my %gs_parpval = ();
my (%gs_idx, $par_sd, $zscore, $pvalue); my $count = 0;
# foreach my $gs (keys %gs_genes) {
foreach my $gsi (@sel_colidx) {
    my $gs = $all_col[$gsi];
    if(($gs_size{$gs} < $imin) or ($gs_size{$gs} > $imax)) { next; }

    $count++; $gs_idx{$gs} = $count;
    @gsref = get_gsscore($gs_genes{$gs}, $gs_size{$gs},
        \%gene_score, \@sel_colidx);

    print ZMAT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    foreach my $col (@sel_colidx) {
        $gs_score{$col}{$gs} = $gsref[0]->{$col};
        $gs_modsize{$col}{$gs} = $gsref[1]->{$col};

        $par_sd = ($sd_gene_score{$col} / sqrt($gs_size{$gs}));
        $zscore = ($gs_score{$col}{$gs} - $mean_gene_score{$col}) / $par_sd;
		print ZMAT "\t", sprintf("%.3f", $zscore);

        # if($iqval) { $gs_zscore{$col}{$gs} = $zscore; }
        $gs_zscore{$col}{$gs} = $zscore;

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

# Multiple Hypotheses Testing correction
$time = runtime(); print "$time: Multiple hypothesis correction ...\n";

my %gs_paresfdr = (); my %gs_sigcount = ();
foreach my $col (@sel_colidx) {
    $gs_paresfdr{$col} = get_esfdr($gs_parpval{$col});

    unless($iqval) { next; }

    foreach my $gs (keys %{$gs_paresfdr{$col}}) {
        if($gs_paresfdr{$col}{$gs} >= (-1)*log($iqval)/log(10)) {
            $gs_sigcount{$gs}++; } }
}

# Printing results
if($iqval) {
    $time = runtime(); print "$time: Printing trimmed Z-score matrix ...\n";

    (my $otzmat = $ozmat) =~ s/\.zscore/\.trim\.zscore/g;
    open TZMAT, ">$otzmat";

    open ZMAT, "$ozmat" or die "Can't open $ozmat!";
    while (<ZMAT>) {
        if($_ =~ /^GS\.ID/) { print TZMAT "$_"; next; }
        chomp($_); @p = split '\t', $_;
        if($gs_sigcount{$p[0]} == 0) { next; }

        print TZMAT "$_\n";
    }

    close ZMAT; close TZMAT;
}
else {
    $time = runtime(); print "$time: Printing ES-FDR matrix ...\n";

    my $sign;
    foreach my $gs (sort {$gs_idx{$a} <=> $gs_idx{$b}} keys %gs_idx) {
        print QMAT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
        foreach my $col (@sel_colidx) {
            $sign = 1;
            unless($gs_score{$col}{$gs} == 0) {
                $sign = ($gs_score{$col}{$gs} / abs($gs_score{$col}{$gs})); }
            print QMAT "\t", $sign*$gs_paresfdr{$col}{$gs};

            if($all_col[$col] eq $gs) {
                print TAB "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
                print TAB "\t", sprintf("%.3f", $gs_zscore{$col}{$gs});
                print TAB "\t", $sign*$gs_paresfdr{$col}{$gs}, "\n"; } }
        print QMAT "\n"; }

    close QMAT;
    close TAB;
}

$time = runtime(); print "$time: DONE ...\n\n";

# Subroutines
# ============
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
            $size{$gs}++;
        }
    }
    close GMT;

    return (\%genes, \%desc, \%size);
}

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
                $mod_size{$col}--;
            }
            else { $totscore += $scoreref->{$col}->{$gene}; }
        }

        $mean_score{$col} = ($totscore/$mod_size{$col});
    }
    
    return (\%mean_score, \%mod_size);
}

# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $pvalues_ref = shift;
    
    my $min_pval = 1; my $num_tests = 0;
    foreach my $gs (keys %{$pvalues_ref}) {
        $num_tests++;
        if(($pvalues_ref->{$gs} != 0) 
                and ($pvalues_ref->{$gs} < $min_pval)) {
            $min_pval = $pvalues_ref->{$gs};
        }
    }

    my %esfdr = (); my $rank = $num_tests;
    foreach my $gs (sort {$pvalues_ref->{$b} <=> $pvalues_ref->{$a}}
        keys %{$pvalues_ref}) {
        if($pvalues_ref->{$gs} == 0) {
            $pvalues_ref->{$gs} = 0.1*sprintf("%.6g", $min_pval);
        }

        $esfdr{$gs} = sprintf("%.3f",
            (-1)*log(($pvalues_ref->{$gs} * $num_tests) / $rank)/log(10));

        if($esfdr{$gs} < 0) { $esfdr{$gs} = 0; }

        $rank--;
    }

    return \%esfdr;
}

# Cleans filename
sub clean_fn {
    my $str = shift;
    $str =~ s/\.\.\///g;
    $str =~ s/^.*\///g;
    $str =~ s/\.\w*$//g;
    return $str;
}

__END__

=head1

Performs Parametric Analysis of Geneset Enrichment (PAGE) based on aggregate of
member-gene scores. This version is specifically for the scenario where the
columns uniquely correspond to genesets and the goal is to assess whether the
right geneset is enriched in the right column.

=head1 USAGE

./gsea_PAGE_self.pl [--imat INPUT_MAT] [--icol COLUMNS_LIST] [--igmt GENESETS_GMT] [--help]

=head1 DESCRIPTION

This script takes a matrix containing 'scores' for all profiled genes (rows)
under any number of different factors (columns), and assesses how unlikely are
the aggregate scores for genes in each geneset in the provided collection.

=head1 ARGUMENTS

=over 12

=item C<--imat>

Input matrix. Contains one header row and the first column contains gene ids.
Can contain NAs.

=item C<--icol>

File containing the list of columns/genesets to consider. This order will be
used to format the z-score and esfdr matrices.

=item C<--igmt>

Collection of genesets in GMT format.

=item C<--imin>

(Optional) Min. no. of genes/geneset. Genesets with less than <imin> genes
will be excluded. Cannot be < 10. Default 10.

=item C<--imax>

(Optional) Max. no. of genes/geneset. Genesets with greater than <imax>
genes will be excluded. Default 200.

=item C<--iabs>

(Optional) Take absolute values of the gene scores & carry-out a one-sided test.

=item C<--itail>

(Optional) 'utail' or 'ltail' for calculating lower or upper tail p-value.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2011 Feb 15

=cut

