#!/usr/bin/perl
use strict;
use warnings;
use PDF;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);

# Creates a gs -> genes hash from a .gmt file, or a gene -> gs hash otherwise
sub assign_genes2gs {
    my $file = shift;
    open FH, "$file" or die "Can't open $file!";
    chomp(my @gset=<FH>); close FH;
    
    my ($gs, %gs_desc, %gs_genes, %genes, @p);
    
    foreach (@gset) {
        @p = split '\t', $_;
        $gs = shift(@p);

        $p[0] =~ s/ \([0-9]*\)$//g;
        $gs_desc{$gs} = $p[0];

        shift(@p);
        foreach my $gene (@p) {
            $gs_genes{$gs}{$gene}++;
            $genes{$gene}++;
        }
    }

    return (\%gs_genes, \%gs_desc, \%genes);
}

# Restrict genes in geneset to specified background
sub bring2_commonbg {
    my $gsref = shift;
    my $generef = shift;
    my %size = (); my $min_size = 10; my $max_size = 250;

    foreach my $gs (keys %{$gsref}) {
        $size{$gs} = 0;
        foreach my $gene (keys %{$gsref->{$gs}}) {
            if(exists $generef->{$gene}) {
                $size{$gs}++;
            }
            else {
                delete $gsref->{$gs}->{$gene};
            }
        }
    }

    foreach my $gs (keys %size) {
        if(($size{$gs} < $min_size) or ($size{$gs} > $max_size)) {
            delete $gsref->{$gs};
        }
    }

    return %size;
}

# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $pvalues_ref = shift;
    
    my $min_pval = 1; my $num_tests = 0;
    foreach my $sets (keys %{$pvalues_ref}) {
        $num_tests++;
        if(($pvalues_ref->{$sets} != 0) 
                and ($pvalues_ref->{$sets} < $min_pval)) {
            $min_pval = $pvalues_ref->{$sets};
        }
    }

    my %esfdr = (); my $rank = $num_tests;
    foreach my $sets (sort {$pvalues_ref->{$b} <=> $pvalues_ref->{$a}}
        keys %{$pvalues_ref}) {
        if($pvalues_ref->{$sets} == 0) {
            $pvalues_ref->{$sets} = 0.1*sprintf("%.6g", $min_pval);
        }

        $esfdr{$sets} =
            (-1)*log(($pvalues_ref->{$sets} * $num_tests) / $rank)/log(10);
        if($esfdr{$sets} < 0) { $esfdr{$sets} = 0; }

        $rank--;
    }

    return %esfdr;
}

my($in_gmt1, $in_gmt2, $out_file) = @ARGV;
open OUT, ">$out_file";

my ($gsgene_ref, $gsdec_ref, $genes_ref, $time);

$time = runtime(); print "\n$time: Assigning genes-to-genesets and trimming them to bg ...";
($gsgene_ref, $gsdec_ref, $genes_ref) = assign_genes2gs($in_gmt1);
my %gs1_genes = %{$gsgene_ref};
my %gs1_desc = %{$gsdec_ref};
my %gs1_allgenes = %{$genes_ref};

($gsgene_ref, $gsdec_ref, $genes_ref) = assign_genes2gs($in_gmt2);
my %gs2_genes = %{$gsgene_ref};
my %gs2_desc = %{$gsdec_ref};
my %gs2_allgenes = %{$genes_ref};

my %gs1_size = bring2_commonbg(\%gs1_genes, \%gs2_allgenes);
my %gs2_size = bring2_commonbg(\%gs2_genes, \%gs1_allgenes);

my %all_genes = %gs1_allgenes;
@all_genes{keys %gs2_allgenes} = values %gs2_allgenes;
my $univ = scalar keys %all_genes;

my @gs1_array = sort keys %gs1_genes;
my @gs2_array = sort keys %gs2_genes;

print "\n\tGS1: $in_gmt1\n\t\tNo. genesets: ", scalar @gs1_array, "\n\t\t";
print "No. genes: ", scalar keys %gs1_allgenes, "\n";
print "\tGS2: $in_gmt2\n\t\tNo. genesets: ", scalar @gs2_array, "\n\t\t";
print "No. genes: ", scalar keys %gs2_allgenes, "\n\n";

print OUT "GS1.ID\tGS1.Size\tGS1.Desc\t";
print OUT "GS2.ID\tGS2.Size\tGS2.Desc\t";
print OUT "No.Common\tOvlp\tJac\tLOR\tESQval\n";

$time = runtime(); print "$time: Starting pair-wise comparisons ...";
my ($tag, $size1, $size2, $min, $common, %tempu, $union, $count);
my ($pval, $best_pval, $best_match, $ncom, $ovlp, $jac, $lor);
my %num_common = (); my %pvalue = (); my %lor = ();
my %overlap = (); my %jaccard = ();
GS1: foreach my $gs1 (@gs1_array) {
    $count++;
    $size1 = $gs1_size{$gs1};
    $best_pval = 1;

    GS2: foreach my $gs2 (@gs2_array) {
        if($gs1 eq $gs2) { next GS2; }

        $size2 = $gs2_size{$gs2};
        $min = $size2; if($size1 < $size2) { $min = $size1; }
 
        $common = 0;
        foreach my $gene (keys %{$gs2_genes{$gs2}}) {
            if(exists $gs1_genes{$gs1}{$gene}) {
                $common++;
            }
        }

        %tempu = %{$gs1_genes{$gs1}};
        @tempu{keys %{$gs2_genes{$gs2}}} = values %{$gs2_genes{$gs2}};
        $union = scalar keys %tempu;

        $pval = hypergeometric_tail($univ, $size1, $size2, $common);
        if($pval < $best_pval) {
            $best_pval = $pval; $best_match = $gs2;
            $ncom = $common;
            $ovlp = sprintf("%.3f", ($ncom/$min));
            $jac = sprintf("%.3f", ($ncom/$union));
            $lor = 0; if($ncom > 0) { $lor = sprintf("%.3f",
                    log(($common/$size1)/($size2/$univ))/log(2)); }
        }
    }

    $tag = $gs1.'__'.$best_match;
    $num_common{$tag} = $ncom;
    $overlap{$tag} = $ovlp;
    $jaccard{$tag} = $jac;
    $lor{$tag} = $lor;
    $pvalue{$tag} = $best_pval;

    if(($count % 1000) == 0) {
        $time = runtime();
        print "\n\t$time: $count ...";
    }
}

$time = runtime(); print "\n\n$time: Performing multiple-hypothesis correction ...";
my %esqvalue = get_esfdr(\%pvalue);

$time = runtime(); print "\n$time: Printing results ...";
my ($gs1, $gs2);
foreach my $tag (sort {$esqvalue{$b} <=> $esqvalue{$a}} keys %esqvalue) {
    # if($esqvalue{$tag} < 2) { last; }
    ($gs1, $gs2) = split '__', $tag;

    print OUT "$gs1\t$gs1_size{$gs1}\t$gs1_desc{$gs1}\t";
    print OUT "$gs2\t$gs2_size{$gs2}\t$gs2_desc{$gs2}\t";
    print OUT "$num_common{$tag}\t$overlap{$tag}\t";
    print OUT "$jaccard{$tag}\t$lor{$tag}\t$esqvalue{$tag}\n";
}

close OUT;

$time = runtime(); print "\n$time: DONE\n\n";

