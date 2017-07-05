#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

#################################################################	
# Future improvements:                                          #
# Calculate overall modularity of the network (Newmann & Girvan #
# [2004] PRE):                                                  #
# 	Q = Sum{s=1_to_m} [ (l_s/L) - (d_s/2*L)^2 ],                #
# where the sum is over the m modules of the partition, l_s is  #
# the no. of links inside module s, L is the total no. of links #
# in the network, and d_s is the total conn of the nodes in     #
# modules.                                                      #
# The first term of the summand in the equation is the fraction #
# of links inside module s; the second term, in contrast,       #
# represents the expected fraction of links in that module, if  #
# links were located at random in the# network (under the only  #
# constraint that the conn sequence coincides with the one of   #
# the original graph).                                          #
#################################################################	

# in_file: Hubber output
my($in_file, $out_file) = @ARGV;

open FH, "$in_file" or die "Can't open $in_file!"; chomp(my @f=<FH>); close FH;
open HH, ">$out_file";

print HH "GS\tSize\t";
print HH "Num.In\tSum.In\tAvg.In\t";
print HH "Num.InOut\tSum.InOut\tAvg.InOut\t";
print HH "Segrn\tAvg.Segrn\n";

my @p;
my ($num_in, $num_inout); #, $num_out);
my ($sum_in, $sum_inout); #, $sum_out);
my ($avg_in, $avg_inout); #, $avg_out);
my ($segrn, $avg_segrn);

shift(@f); shift(@f);
foreach (@f) {
    @p = split '\t', $_;

    $num_in = $p[7];
    $avg_in = $p[5];
    $sum_in = ($avg_in*$num_in);

    $num_inout = ($p[4] - $num_in);
    $sum_inout = (($p[4]*$p[2]) - $sum_in);
    if($num_inout != 0) { $avg_inout = ($sum_inout/$num_inout); }
    else { $avg_inout = 'NA'; }

    # $num_out = ($num_inout - $num_in);
    # $sum_out = ($sum_inout - $sum_in);
    # if($num_out != 0) { $avg_out = ($sum_out/$num_out); }
    # else { $avg_out = 'NA'; }

    if($sum_inout != 0) { $segrn = ($sum_in/$sum_inout); }
    else { $segrn = 'NA'; }

    if($avg_inout != 0) { $avg_segrn = ($avg_in/$avg_inout); }
    else { $avg_segrn = 'NA'; }

    print HH "$p[0]\t$p[1]";
    print HH "\t$num_in\t", sprintf("%.3f", $sum_in), "\t", sprintf("%.3f", $avg_in);
    print HH "\t$num_inout\t", sprintf("%.3f", $sum_inout), "\t", sprintf("%.3f", $avg_inout);
    print HH "\t", sprintf("%.3f", $segrn), "\t", sprintf("%.3f", $avg_segrn), "\n";
}

close HH;

