#!/usr/bin/perl
use strict;
use warnings;
use PDF;
#use Getopt::Long;
#use Pod::Usage;
use Data::Dumper;
#use Time::SoFar qw(runtime);

# Calculates median
sub get_median {
    my $hashref = shift;
    my @array = sort {$a <=> $b} grep {!/^NA$/} (values %{$hashref});
    if(@array % 2) { return $array[@array/2]; }
    else { return (($array[(@array/2)-1] + $array[@array/2]) / 2); }
}

my($in_mat, $in_par, $out_colstat) = @ARGV;
open MH, "$in_mat" or die "Can't open $in_mat!"; chomp(my @mat=<MH>); close MH;
open SH, ">$out_colstat";

# common, topovlp, jaccard, log-odds, es-hypergeometric

my @col_names = split '\t', $mat[0]; shift(@col_names);
my $univ = $#mat;

my %col_rows = (); my @p;
for(my $i=1; $i<=$#mat; $i++) {
    @p = split '\t', $mat[$i];
    for(my $j=1; $j<=$#p; $j++) {
        if($p[$j] != 0) {
            $col_rows{$col_names[$j-1]}{$p[0]} = $p[$j];
        }
    }
}

my ($size1, $size2, $min);
my ($pair, $common, %unionhash, $union, $pval);
my (%colpair_ovlp, %colpair_jac, %colpair_lod, %colpair_hyperg);

for(my $j=0; $j<$#col_names; $j++) {
    $size1 = scalar(keys %{$col_rows{$col_names[$j]}});

    for(my $k=($j+1); $k<=$#col_names; $k++) {
        $size2 = scalar(keys %{$col_rows{$col_names[$k]}});
        $min = $size1; if($size2 < $min) { $min = $size2; }

        $pair = join '__', sort($col_names[$j], $col_names[$k]);
        
        $common = 0; $union = 0;
        foreach my $r (keys %{$col_rows{$col_names[$j]}}) {
            if(exists $col_rows{$col_names[$k]}{$r}) {
                $common++;
            }
        }

        %unionhash = ();
        @unionhash{keys %{$col_rows{$col_names[$j]}}}
            = values %{$col_rows{$col_names[$j]}};
        @unionhash{keys %{$col_rows{$col_names[$k]}}}
            = values %{$col_rows{$col_names[$k]}};
        $union = scalar(keys %union);

        $colpair_ovlp{$pair} = $common/$min;
        $colpair_jac{$pair} = $common/$union;
        $colpair_lod{$pair} = ($common/$size1)/($size2/$univ);
        $pval = hypergeometric_tail($univ, $size1, $size2, $common);
        $colpair_hyperg{$pair} = (-1)*log($pval)/log(10);
    }
}

# Std.Dev, IQR, Median Absolute Deviation ...
# http://en.wikipedia.org/wiki/Robust_measures_of_scale S_n, Q_n

my (%median_dist, $median);
for(my $j=0; $j<=$#col_names; $j++) {
    %median_dist = (); $median = 0;
    for(my $k=0; $k<=$#col_names; $k++) {
        if($k == $j) { $next; }
        $pair = join '__', sort($col_names[$j], $col_names[$k]);
        %median_dist{$pair} = $colpair_lod{$pair};
    }
    $median = get_median(\%median_dist);
}

close SH;

