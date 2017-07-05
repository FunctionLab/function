#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use IO::Tee;
use Statistics::Distributions;
use List::Util 'shuffle';
use Time::SoFar qw(runtime);

# Creates array of indices from a string representation
sub get_idxarray {
    my $string = shift;
    my $first_line = shift;
    
    my @ids_array = split("\t", $first_line);
    my @idx_array = (); my (@p, @q);
    
    if($string eq 'all') {
        for(my $j=1;$j<=$#ids_array;$j++) {
            push(@idx_array, $j);
        }
    }
    else {
        @p = split '_', $string;

        foreach(@p) {
            if($_ =~ /:/) {
                @q = split ':', $_;
                for(my $j=$q[0]; $j<=$q[1]; $j++) {
                    push(@idx_array, $j);
                }
            }
            else {
                push(@idx_array, $_);
            }
        }
    }

    return (\@ids_array, \@idx_array);
}

# Calculates median
sub get_median {
    my $hashref = shift;
    my @array = sort {$a <=> $b} grep {!/^NA$/} (values %{$hashref});
    if(@array % 2) { return $array[@array/2]; }
    else { return (($array[(@array/2)-1] + $array[@array/2]) / 2); }
}

# Creates a gs -> genes hash from a .gmt file
sub assign_genes2gs {
    my $file = shift;
    my $bgref = shift;
    
    my (%gs2desc, %gs2genes, %gs2size, @p, $gs);

    print "\n\t*$file*";
    open FH, "$file" or die "Can't open $file!";
        chomp(my @gset=<FH>); close FH;

    foreach my $line (@gset) {
        @p = split '\t', $line;
        $gs = shift(@p);
        $p[0] =~ s/ \([0-9]*\)$//g;
        $gs2desc{$gs} = $p[0];
        shift(@p);

        foreach my $g (@p) {
            unless(exists $bgref->{$g}) { next; }
            push(@{$gs2genes{$gs}}, $g);
            $gs2size{$gs}++;
        }
    }

    return (\%gs2genes, \%gs2desc, \%gs2size);
}

# Calclates mean gs score in all chosen columns
sub get_gsscore {
    my $gs_generef = shift;
    my $size = shift;
    my $gene_scoreref = shift;
    my $colref = shift;
    
    my %mean_score=(); my %mod_size=(); my $col;

    for(my $k=0;$k<=$#{$colref};$k++) {
        $col = ${$colref}[$k];

        my $totscore = 0;
        my $mod_size{$col} = $size;

        foreach my $gene (@{$gs_generef}) {
            if($gene_scoreref->{$col}->{$gene} eq 'NA') { $mod_size{$col}--; }
            else { $totscore += $gene_scoreref->{$col}->{$gene}; }
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

    my $fdr; my %esfdr = (); my $rank = $num_tests;
    foreach my $gs (sort {$pvalues_ref->{$b} <=> $pvalues_ref->{$a}}
        keys %{$pvalues_ref}) {
        if($pvalues_ref->{$gs} == 0) {
            $pvalues_ref->{$gs} = 0.1*sprintf("%.6g", $min_pval);
        }

        $fdr = (($pvalues_ref->{$gs} * $num_tests) / $rank);
        if($fdr > 1) { $fdr = 1; }

        $esfdr{$gs} = sprintf("%.3f", (-1)*log($fdr)/log(10));

        $rank--;
    }

    return \%esfdr;
}

# Reading in files and options
$min_genes = 10; $max_genes = 200;
my($help, $time, $in_mat, $in_coltag, $in_gset);

my @col_array = get_idxarray($coltag, $mat[0]); # One more array is returned by this function

# Populating genes scores and background population
my %num_genes = (); my %mean_gene_score = (); my %sumsq_gene_score = ();
my $col;
for(my $j=0; $j<=$#col_array; $j++) {
    $col = $col_array[$j];
    $num_genes{$col} = 0;
    $mean_gene_score{$col} = 0;
    $sumsq_gene_score{$col} = 0;
}

# Mean and SumSq calculation using Welford's algorithm
my %gene_score = (); my %bg_genes = (); my @p;
my ($score, $old_mean, $old_var);
for(my $i=1; $i<=$#mat; $i++) {
    if($mat[$i] =~ /^#/) { next; }

    $mat[$i] =~ s///g;
    @p=split("\t",$mat[$i]);

    $bg_genes{$p[0]}++;
    for(my $j=0; $j<=$#col_array; $j++) {
        $col = $col_array[$j]; $score = $p[$col_array[$j]];

        $gene_score{$col}{$p[0]} = $score;

        unless($score eq 'NA') {
            $num_genes{$col}++;
            
            $old_mean = $mean_gene_score{$col};
            $mean_gene_score{$col} = $old_mean +
                ($score-$old_mean)/$num_genes{$col};

            $old_var = $sumsq_gene_score{$col};
            $sumsq_gene_score{$col} = $old_var +
                ($score-$old_mean)*($score-$mean_gene_score{$col});
        }
    }
}

# Calculating sd and median
my %sd_gene_score = (); my %median_gene_score = ();
for(my $j=0; $j<=$#col_array; $j++) {
    $col = $col_array[$j];
    $sd_gene_score{$col} = sqrt($sumsq_gene_score{$col}/$num_genes{$col});
    $median_gene_score{$col} = get_median($gene_score{$col});
}

# Populating genests
my @gsref = assign_genes2gs{$in_gset, \%bg_genes};
my %gs_genes = %{$gsref[0]};
my %gs_desc = %{$gsref[1]};
my %gs_size = %{$gsref[2]};

# Recording geneset scores, and printing zscores
my %gs_score = (); my %gs_modsize = (); my %gs_parpval = ();
my %gs_idx, @gs_scoreref, $par_sd, $zscore, $pvalue; my $count = 0;
foreach my $gs (keys %gs_genes) {
    if(($gs_size{$gs} < $min_genes) or ($gs_size{$gs} > $max_genes)) { next; }

    $count++; $gs_idx{$gs} = $count;
    @gs_scoreref = get_gsscore($gs_genes{$gs}, $gs_size{$gs},
        \%gene_score, \@col_array);

    print HH "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    for(my $j=0; $j<=$#col_array; $j++) {
        $col = $col_array[$j];
        $gs_score{$col}{$gs} = $gs_scoreref[0]->{$col};
        $gs_modsize{$col}{$gs} = $gs_modsizeref[0]->{$col};

        $par_sd = ($sd_gene_score{$col}/sqrt($gs_size{$gs}));
        $zscore = sprintf("%.3f",
            ($gs_score{$col}{$gs}-$mean_gene_score{$col})/$par_sd);
        print HH "\t$zscore";

        $pvalue = sprintf("%.5g",
            (Statistics::Distributions::uprob(abs($zscore)))*2);
        $gs_parpval{$col}{$gs} = $pvalue;

        # $old_mean
    }
    print HH "\n";
}

# Creating random background gs scores
my $rand_trials = 1000; my @rand_genearray;
my %rand_gsscore = (); my %emp_mean = (); my %emp_sd = ();
my ($old_mean, $curr_mean); # my $old_emp_mean = 0; my $old_emp_sd = 0;
for(my $k=1; $k<=$rand_trials; $k++) {
    @rand_genearray = shuffle(keys %bg_genes);

    for(my $j=0; $j<=$#col_array; $j++) {
        $col = $col_array[$j]; $s=3;
        
        $old_mean =
            ($rand_genearray[0]+$rand_genearray[1]+$rand_genearray[2])/3;
        push(@{$rand_gsscore{$col}{$s}}, $old_mean);
        $old_emp_mean = $old_mean;

        for(my $s=4; $s<10; $s++) {
            $curr_mean = $old_mean + ($rand_genearray[$s-1]-$old_mean)/$s;

            push(@{$rand_gsscore{$col}{$s}}, $curr_mean);
            $emp_mean{$col}{$s} = $old_emp_mean + ($curr_mean-$old_emp_mean)/$k;
        }
    }
}

# Multiple hypothesis correction
my %gs_paresfdr = ();
for(my $j=0; $i<=$#col_array; $j++) {
    $col = $col_array[$j];
    $gs_paresfdr{$col} = get_esfdr($gs_parpval{$col});
}

# Printing FDR
my $sign;
foreach my $gs (sort {$gs_idx{$a} <=> $gs_idx{$b}} keys %gs_paresfdr) {
    print JH "$gs\t$gs_desc{$gs}\t$gs_size{$gs}";
    for(my $j=0; $j<=$#col_array; $j++) {
        $col = $col_array[$j];
        $sign = ($gs_score{$col}{$gs}/abs($gs_score{$col}{$gs}));
        print JH "\t", $sign*$gs_paresfdr{$col}{$gs};
    }
    print JH "\n";
}



