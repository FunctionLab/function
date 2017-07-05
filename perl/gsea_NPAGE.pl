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


my ($help, $imat, $icol, $igmt, $ibgc, $osmat, $oqmat, $iabs, $iqval);
my $imin = 3; my $imax = 200; my $intv = 1000; my $irnd = 10000;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'icol=s' => \$icol,
          'ibgc=s' => \$ibgc,
          'intv=i' => \$intv,
          'igmt=s' => \$igmt,
          'imin=i' => \$imin,
            'iabs' => \$iabs,
          'imax=i' => \$imax,
          'irnd=i' => \$irnd,
         'iqval=f' => \$iqval) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);
if($imin < 3) { $imin = 3; print "\nimin set to 3.\n"; }
if($imin > 200) { $imax = 200; print "\nimax set to 200.\n"; }

# Reading in options and all files
my $time = runtime();
print "\n$time: Reading in options and files ...";

open FH,"$imat"; chomp(my @mat=<FH>); close FH;

my $otag = ''; if($iabs) { $otag = '.abs'; }
$osmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.score.mat';
$oqmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.esfdr.mat';
open SMAT,">$osmat"; open QMAT,">$oqmat";

unless($icol) { $icol = 'all'; }
my ($idsref, $idxref) = get_idxarray($icol, shift(@mat));
my @all_col = @{$idsref}; my @sel_colidx = @{$idxref};

my $numcol = scalar(@sel_colidx);
print "\nInput file: $imat\n";
print "No. rows: $#mat\tNo. cols: $numcol\n";

# Printing column names and indexing them
$time = runtime();
print "\n$time: Indexing input genes, genesets desc and background lists ...\n";

print QMAT "GS.ID\tGS.Desc\tGS.Size";
print SMAT "GS.ID\tGS.Desc\tGS.Size";

my %num_genes = ();
foreach my $col (@sel_colidx) {
    $num_genes{$col} = 0;

    print QMAT "\t$all_col[$col]";
    print SMAT "\t$all_col[$col]"; }

print QMAT "\n";
print SMAT "\n";


# Indexing background genes and their associated scores
my %bg_genes = my %gene_intv = my %intv_genes = ();
if($ibgc) {
    open BG, "$ibgc" or die "Can't open $ibgc!";
    while (<BG>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $bg_genes{$p[0]} = $p[1]; }
    close BG;

    my @genes = sort {$bg_genes{$a} <=> $bg_genes{$b}} keys %bg_genes;
    my $N = scalar @genes; my $k = 1;
    my @quant = (); my $q;

    for(my $i=0; $i<=$#genes; $i++) {
        if($bg_genes{$genes[$i]} <= ($N*$k/$intv)) {
            $gene_intv{$genes[$i]} = $k;
            $intv_genes{$k}{$genes[$i]}++; }
        else { $i--; $k++; } } }


# Indexing gene scores and no. genes per column
my %gene_score = ();
foreach (@mat) {
    if($_ =~ /^#/) { next; }
    my @p = split '\t', $_;
    if($ibgc) { unless(exists $bg_genes{$p[0]}) { next; } }
    else { $bg_genes{$p[0]}++; }

    foreach my $col (@sel_colidx) {
        if($iabs) {
            if(looks_like_number($p[$col])) {
                $p[$col] = abs($p[$col]); } }
        $gene_score{$col}{$p[0]} = $p[$col];

        if(looks_like_number($p[$col])) {
            $num_genes{$col}++; } } }


# Indexing gs descriptions and member genes
$time = runtime(); print "$time: Indexing gs descriptions and genes ...\n";

my @gsref = parse_gmt($igmt, \%bg_genes);
my %gs_genes = %{$gsref[0]};
my %gs_desc = %{$gsref[1]};
my %gs_size = %{$gsref[2]};


# Recording geneset scores
$time = runtime();
print "$time: Calculating gs scores ...", "\n";

my %gs_score = my %gs_modsize = my %gs_idx = my %col_sizegs = ();
my $count = 0;
foreach my $gs (sort keys %gs_genes) {
    if(($gs_size{$gs} < $imin) or ($gs_size{$gs} > $imax)) { next; }

    $count++; $gs_idx{$gs} = $count;
    my ($score, $size) = get_gsscore($gs_genes{$gs}, \%gene_score,
        \@sel_colidx);

    foreach my $col (@sel_colidx) {
        $gs_score{$col}{$gs} = $score->{$col};
        $gs_modsize{$col}{$gs} = $size->{$col};

        $col_sizegs{$col}{$size->{$col}}{$gs}++; } }


# Calculating empirical significance scores
my %gs_pval = ();
foreach my $col (@sel_colidx) {
    foreach my $size (sort keys %{$col_sizegs{$col}}) {
        my %pval;
        if($ibgc) {
            %pval = emp_dist($gene_score{$col}, $gs_score{$col}, $size,
                $col_sizegs{$col}{$size}, \%gene_intv, \%intv_genes); }
        else {
            %pval = emp_dist($gene_score{$col}, $gs_score{$col}, $size,
                $col_sizegs{$col}{$size}); }

        foreach my $gs (keys %pval) {
            $gs_pval{$col}{$gs} = $pval{$gs}; } } }


# Multiple Hypotheses Testing correction
$time = runtime(); print "$time: Multiple hypothesis correction ...\n";

my %gs_esfdr = ();
foreach my $col (@sel_colidx) {
    $gs_esfdr{$col} = get_esfdr($gs_pval{$col}); }


# Printing results
$time = runtime(); print "$time: Printing results ...\n";

my $sign;
foreach my $gs (sort {$gs_idx{$a} <=> $gs_idx{$b}} keys %gs_idx) {
    print SMAT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    print QMAT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    foreach my $col (@sel_colidx) {
        print SMAT "\t$gs_score{$col}{$gs}";

        $sign = ($gs_score{$col}{$gs} / abs($gs_score{$col}{$gs}));
        print QMAT "\t", $sign*$gs_esfdr{$col}{$gs}; }
    print SMAT "\n";
    print QMAT "\n"; }

$time = runtime(); print "$time: DONE ...\n\n";

close SMAT;
close QMAT;


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
        $desc{$gs} = shift(@p);

        foreach my $g (@p) {
            unless(exists $bgref->{$g}) { next; }
            $genes{$gs}{$g}++;
            $size{$gs}++; } }
    close GMT;

    return (\%genes, \%desc, \%size); }

# Calclates mean gs score in all chosen columns
sub get_gsscore {
    my $generef = shift;
    my $scoreref = shift;
    my $colref = shift;
    
    # my $size = scalar keys %{$generef};
    my %mean_score = my %mod_size = ();

    foreach my $col (@{$colref}) {
        my $tsize = my $tmean = 0;

        foreach my $g (keys %{$generef}) {
            my $s = $scoreref->{$col}->{$g};
            unless(looks_like_number($s)) { next; }
            $tsize++;
            $tmean += ($s - $tmean) / $tsize; }

        $mod_size{$col} = $tsize;
        $mean_score{$col} = $tmean; }

    return (\%mean_score, \%mod_size); }

# Calculates empirical distribution and p-value
sub emp_dist {
    my $gene_scoreref = shift; # hash of gene scores
    my $gs_scoreref = shift; # hash of gs scores
    my $size = shift; # geneset size
    my $gsref = shift; # hash of genesets

    my $intg = shift;
    my $gint = shift;

    my %remain_gs = %{$gsref};

    # my @permut_values = ();
    my (@randg, $rmean, $pval, $ci, %exc, %pvalue);

    foreach my $gs (keys %remain_gs) {
        $exc{$gs} = 1; $pvalue{$gs} = 1/$irnd; }

    for(my $r=1; $r<$irnd; $r++) {
        if($ibgc) {
            @randg = ();
            foreach my $g (keys %{$gene_scoreref}) {
                my $int = $gint->{$g};
                my @rgenes = shuffle grep {!/^$g$/} keys %{$intg->$int};
                push(@randg, $rgenes[0]); } }
        else {
            @randg = (shuffle keys %{$gene_scoreref})[0..($size-1)]; }
        
        $rmean = 0;
        for(my $s=0; $s<$size; $s++) {
            $rmean += (($gene_scoreref->{$randg[$s]} - $rmean) / ($s + 1)); }

        # push(@permut_values, $rmean);

        foreach my $gs (keys %remain_gs) {
            if(abs($rmean) >= abs($gs_scoreref->{$gs})) { $exc{$gs}++; }
            if($r < 10) { next; }

            $pval = ($exc{$gs} / $r);

            if($pval >= (10 / $r)) {
                $pvalue{$gs} = $pval;
                delete $remain_gs{$gs}; } }

        if(scalar keys %remain_gs == 0) { last; } }

    # return (\@permut_values, \%pvalue);
    return \%pvalue; }

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

        $rank--; }

    return \%esfdr; }

# Cleans filename
sub clean_fn {
    my $str = shift;
    $str =~ s/~//g;
    $str =~ s/\.\.\///g;
    $str =~ s/^.*\///g;
    $str =~ s/\.\w*$//g;
    return $str; }

__END__

=head1

Performs Non-Parametric Analysis of Geneset Enrichment (NPAGE) based on aggregate of
member-gene scores.

=head1 USAGE

./gsea_NPAGE.pl [--imat INPUT_MAT] [--icol COLUMNS] [--igmt GENESETS_GMT] [--help]

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

(Optional) String to specify columns of interest using ':' to specific a continuous series
& '_' to specify a break. Index begins at 0. E.g. "3_5:8_11_13" would mean
columns 3, 5, 6, 7, 8, 11 & 13. Default is to use all columns.

=item C<--igmt>

Collection of genesets in GMT format.

=item C<--imin>

(Optional) Min. no. of genes/geneset. Genesets with less than <imin> genes
will be excluded. Cannot be < 3. Default 3.

=item C<--imax>

(Optional) Max. no. of genes/geneset. Genesets with greater than <imax>
genes will be excluded. Default 200.

=item C<--imax>

(Optional) No. of iterations of permutation to run. Default 10000.

=item C<--iabs>

(Optional) If this option is provided, absolute values of the gene scores are
considered & a one-sided test is carried out.

=item C<--iqval>

(Optional) If this option is provided, then only the Z-score matrix is printed,
with the rows (genesets) with a Qvalue <= iqval in at least one column.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2011 Feb 15

=cut

