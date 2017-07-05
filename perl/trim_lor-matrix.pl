#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my($in_mat, $in_tab, $out_mat) = @ARGV;
open MH, "$in_mat" or die "Can't open $in_mat!";
    chomp(my @mat=<MH>); close MH;
open TH, "$in_tab" or die "Can't open $in_tab!";
    chomp(my @tab=<TH>); close TH;
open HH, ">$out_mat";

my %sel_rows = (); my %sel_cols = (); my %sel_pairs = ();
shift(@tab); my (@p, $tag);
foreach (@tab) {
    @p = split '\t', $_;
    $sel_rows{$p[0]} = $p[2].' ('.$p[1].')';
    $sel_cols{$p[3]} = $p[5].' ('.$p[4].')';

    $tag = $p[0].'__'.$p[3];
    $sel_pairs{$tag}++;
}

my @row_array = ();
my @col_array = split '\t', shift(@mat); shift(@col_array); shift(@col_array);

my $tot_rows = scalar (@mat) - 2; my $tot_cols = scalar(@col_array);
my $num_sel_rows = scalar keys %sel_rows;
my $num_sel_cols = scalar keys %sel_cols;
print "\nNo. total rows: $tot_rows\nNo. total cols: $tot_cols\n";
print "\nNo. select. rows: $num_sel_rows\nNo. select. cols: $num_sel_cols\n\n";

shift(@mat);

my $i = 0; my @lor = ();
my %row_num_nonzero = (); my %col_num_nonzero = ();

foreach (@mat) {
    @p = split '\t', $_;
    push(@row_array, shift(@p)); shift(@p);
    
    for(my $j=0; $j<=$#p; $j++) {
        $tag = $row_array[$#row_array].'__'.$col_array[$j];
        if(exists $sel_pairs{$tag}) {
            $lor[$i][$j] = $p[$j];
            $row_num_nonzero{$i}++;
            $col_num_nonzero{$j}++;
        }
        else { $lor[$i][$j] = 0 }
    }
    $i++;
}

print HH "\t";
for(my $j=0; $j<=$#col_array; $j++) {
    unless(exists $col_num_nonzero{$j}) { next; }
    print HH "\t$col_array[$j]";
}
print HH "\n\t";
for(my $j=0; $j<=$#col_array; $j++) {
    unless(exists $col_num_nonzero{$j}) { next; }
    print HH "\t$sel_cols{$col_array[$j]}";
}
print HH "\n";

ROW: for(my $i=0; $i<=$#row_array; $i++) {
    unless(exists $row_num_nonzero{$i}) { next ROW; }
    print HH "$row_array[$i]\t$sel_rows{$row_array[$i]}";

    COL: for(my $j=0; $j<=$#col_array; $j++) {
        unless(exists $col_num_nonzero{$j}) { next COL; }
        print HH "\t$lor[$i][$j]";
    }
    print HH "\n";
}

close HH;
