#!/usr/bin/perl

# ================================================
# Name : 
# Purpose : 
# Created : 27-01-2015
# Last Modified : Tue 27 Jan 2015 03:17:36 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

use String::Approx qw(adist adistr);
print "\n", (1 - (abs adistr($ARGV[0], $ARGV[1]))), "\n";

