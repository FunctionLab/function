#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

# This script currently takes in a machine-learning dataset in the format below
# and replaces every missing value (i,j) with the median of the (,j) values
# calculated only using instances from the same class.
# Format (each line): <label> <index1>:<value1> <index2>:<value2> ...

# Calculates median
sub get_median {
    my $hashref = shift;
    my @array = sort {$a <=> $b} grep {!/^NA$/} (values %{$hashref});
    if(@array % 2) { return $array[@array/2]; }
    else { return (($array[(@array/2)-1] + $array[@array/2]) / 2); }
}

my $in_file = $ARGV[0];
open FH, "$in_file" or die "Can't open $in_file!"; chomp(my @f=<FH>); close FH;
my $out_file = $in_file.'.imputed';
open HH, ">$out_file";

my $i = 0; my ($feature, $class, $fc, $value, @p); my $max = -1000; my $min = 1000;
my %all_instances = (); my %all_features = (); my %fc_values = ();
foreach (@f) {
    $i++;
    @p = split ' ', $_;
    $class = shift(@p);
    $all_instances{$i} = $class;

    foreach my $f (@p) {
        ($feature, $value) = split ':', $f;
        $all_features{$feature}++;
        if($value eq 'NA') { next; }
        $fc = $feature.'__'.$class;
        $fc_values{$fc}{$i} = $value;

        if($value > $max) { $max = $value; }
        if($value < $min) { $min = $value; }
    }
}

my %fc_median = ();
foreach my $fc (keys %fc_values) {
    $fc_median{$fc} = get_median($fc_values{$fc});
}

foreach my $i (sort {$a <=> $b} keys %all_instances) {
    $class = $all_instances{$i};
    print HH "$class";
    
    foreach my $feature (sort {$a <=> $b} keys %all_features) {
        print HH " $feature:";
        $fc = $feature.'__'.$class;
        if(exists $fc_values{$fc}{$i}) { print HH "$fc_values{$fc}{$i}"; }
        else { print HH "$fc_median{$fc}"; }
    }
    print HH "\n";
}

close HH;

print "\nMax: $max\tMin: $min\n\n";

