#!/usr/bin/perl

# ================================================
# Name : hypergeometric_pvalue.pl
# Purpose : Get hypergeometric pvalue
# Created : 19-09-2013
# Last Modified : Fri 14 Feb 2014 04:56:26 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use lib '/home/arjunk/software/lib/perl5/';
use PDF;

#$univ = 1973; $size1 = 11; $size2 = 5; $ncom = 2;
my ($univ, $size1, $size2, $ncom) = @ARGV;
print hypergeometric_tail($univ, $size1, $size2, $ncom), "\n";

