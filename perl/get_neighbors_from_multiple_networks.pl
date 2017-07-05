#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);

my($in_genelist, $in_netlist, $out_mat) = @ARGV;
open GH, "$in_genelist" or die "Can't open $in_genelist!";
    chomp(my @f=<GH>); close GH;
open NH, "$in_netlist" or die "Can't open $in_netlist!";
    chomp(my @g=<NH>); close NH;
open MH, ">$out_file";

# ...

close MH;
