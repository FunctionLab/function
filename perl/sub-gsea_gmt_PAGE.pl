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

my ($help, $imat, $icol, $igmt, $ozmat, $oqmat, $iabs); my $imin = 10; my $imax = 200;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'icol=s' => \$icol,
          'igmt=s' => \$igmt,
          'imin=i' => \$imin,
            'iabs' => \$iabs,
          'imax=i' => \$imax) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);
if($imin < 10) { print "\nMinimum geneset size cannot be < 10.\n\n"; exit; }

# Reading in options and all files
my $time = runtime();
print "\n$time: Reading in options and files ...";

open FH,"$imat"; chomp(my @mat=<FH>); close FH;

my $otag = ''; if($iabs) { $otag = '.abs'; }
$ozmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.zscore.mat';
$oqmat = clean_fn($imat).$otag.'-'.clean_fn($igmt).'.esfdr.mat';
open ZMAT,">$ozmat"; open QMAT,">$oqmat";

unless($icol) { $icol = 'all'; }
my ($idsref, $idxref) = get_idxarray($icol, shift(@mat));
my @all_col = @{$idsref}; my @sel_colidx = @{$idxref};

my $numcol = scalar(@sel_colidx);
print "\nInput file: $imat\n";
print "No. rows: $#mat\tNo. cols: $numcol\n";

# Populating genes scores and background population
$time = runtime();
print "\n$time: Indexing input genes, genesets desc and background lists ...\n";

print QMAT "GS.ID\tGS.Desc\tGS.Size";
print ZMAT "GS.ID\tGS.Desc\tGS.Size";

my %num_genes = (); my %mean_gene_score = (); my %sumsq_gene_score = ();
foreach my $col (@sel_colidx) {
    $num_genes{$col} = 0;
    $mean_gene_score{$col} = 0;
    $sumsq_gene_score{$col} = 0;

    print QMAT "\t$all_col[$col]";
    print ZMAT "\t$all_col[$col]";
}

print QMAT "\n";
print ZMAT "\n";

# Mean and SumSq calculation using to Welford's algorithm
my %gene_score = (); my %bg_genes = (); my (@p, $old_mean);
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
                ($p[$col] - $old_mean)*($p[$col] - $mean_gene_score{$col});
        }
    }
}

# Calculating sd
my %sd_gene_score = ();
foreach my $col (@sel_colidx) {
    $sd_gene_score{$col} = sqrt($sumsq_gene_score{$col}/$num_genes{$col});
}

# Indexing gs descriptions and member genes
$time = runtime(); print "$time: Indexing gs descriptions and genes ...\n";

my @gsref = parse_gmt($igmt, \%bg_genes);
my %gs_genes = %{$gsref[0]};
my %gs_desc = %{$gsref[1]};
my %gs_size = %{$gsref[2]};

# Recording geneset scores, and printing zscores
$time = runtime();
print "$time: Calculating gs scores ...", "\n";

my %gs_score = (); my %gs_modsize = (); my %gs_parpval = ();
my (%gs_idx, $par_sd, $zscore, $pvalue); my $count = 0;
foreach my $gs (keys %gs_genes) {
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

        if($iabs) { $pvalue = sprintf("%.5g",
                (Statistics::Distributions::uprob($zscore))); }
        else { $pvalue = sprintf("%.5g",
                (Statistics::Distributions::uprob(abs($zscore)))*2); }
        $gs_parpval{$col}{$gs} = $pvalue;
    }
    print ZMAT "\n";
}

# Multiple Hypotheses Testing correction
$time = runtime(); print "$time: Multiple hypothesis correction ...\n";

my %gs_paresfdr = ();
foreach my $col (@sel_colidx) {
    $gs_paresfdr{$col} = get_esfdr($gs_parpval{$col});
}

# Printing results
$time = runtime(); print "$time: Printing results ...\n";

my $sign;
foreach my $gs (sort {$gs_idx{$a} <=> $gs_idx{$b}} keys %gs_idx) {
    print QMAT "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    foreach my $col (@sel_colidx) {
        $sign = ($gs_score{$col}{$gs} / abs($gs_score{$col}{$gs}));
        print QMAT "\t", $sign*$gs_paresfdr{$col}{$gs};
    }
    print QMAT "\n";
}

$time = runtime(); print "$time: DONE ...\n\n";

close QMAT;
close ZMAT;

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
            push(@sel_idx, $j);
        }
    }
    else {
        @p = split '_', $string;
        foreach(@p) {
            @q = split ':', $_;
            push(@sel_idx, $q[0]);
            if($#q eq 0) { next; }

            while($q[0]<$q[1]) {
                $q[0]++; push(@sel_idx, $q[0]);
            }
       }
    }

    return (\@all_ids, \@sel_idx);
}

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
member-gene scores of subset of genesets.

=head1 USAGE

sub-gsea_gmt_PAGE.pl [--imat INPUT_MAT] [--icol COLUMNS] [--igmt GENESETS_GMT] [--help]

=head1 DESCRIPTION

This script takes a matrix containing 'scores' for all profiled genes (rows)
under any number of different factors (columns), and assesses how unlikely are
the aggregate scores for subsets of genes in each geneset in the provided collection.

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
will be excluded. Cannot be less than 10. Default 10.

=item C<--imax>

(Optional) Max. no. of genes/geneset. Genesets with greater than <imax>
genes will be excluded. Default 200.

=item C<--iabs>

(Optional) If this option is provided, absolute values of the gene scores are
considered & a one-sided test is carried out.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2011 Feb 15

=cut

