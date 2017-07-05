#!/usr/bin/perl

# ================================================
# Name : test.pl
# Purpose : 
# Created : 24-12-2014
# Last Modified : Wed 24 Dec 2014 02:35:36 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

chomp(my $l = `wc -l test.dat`); $l =~ s/ .*$//g;
my $sl = int($l/5);

for(my $i=0; $i<100; $i++) {
    `perl -MList::Util -e 'print List::Util::shuffle <>' test.dat | head -n $sl >> test.samp.dat`;
    #`Dat2Dab -i test.dat -u 0.2 >> test.samp.dat`; }
}

