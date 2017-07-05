#!/usr/bin/perl

use strict;
use warnings;

my ($itxt, $otxt) = @ARGV;

open OT, ">$otxt";
open IT, "$itxt" or die "Can't open $itxt!";
while (<IT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    # my @q = split ' ', $p[1];
    #my $gs = shift @p;
    #my $d = shift @p;
    foreach my $w (@p[2..$#p]) {
        # print OT "$gs","_$d\t$g\n"; } }
        # print OT "$gs\t$g\n"; } }
        print OT "$p[0]\t$p[1]\t$w\n"; } }
close IT;
close OT;
