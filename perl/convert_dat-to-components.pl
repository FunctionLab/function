#!/usr/bin/perl
use strict;
use warnings;
use Graph::Undirected;
use List::Util 'shuffle';

my ($idat, $icut, $occs) = @ARGV;


my $graph = Graph::Undirected->new;
my %edges = ();

open DAT, "$idat" or die "Can't open $idat!";
while (<DAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    if(($#p eq 0) or ($p[1] eq '')) {
        $graph->add_vertex($p[0]); }
    else {
        my $e = join '__', sort($p[0], $p[1]);
        $edges{$e} = $p[2];
        if($p[2] >= $icut) {
            $graph->add_edge($p[0], $p[1]); }
        else {
            $graph->add_vertex($p[0]);
            $graph->add_vertex($p[1]); } } }
close DAT;


my @cc = $graph->connected_components();

open CC, ">$occs";
my $count = 0;
foreach my $c (@cc) {
    $count++;
    print CC "cc:", sprintf("%03d", $count);
    print CC "\t", scalar @{$c};
    #@{$c} = shuffle @{$c};

    my %deg = ();
    for(my $i=0; $i<$#{$c}; $i++) {
        for(my $j=($i+1); $j<=$#{$c}; $j++) {
            my $e = join '__', sort(${$c}[$i], ${$c}[$j]);
            $deg{${$c}[$i]} += $edges{$e};
            $deg{${$c}[$j]} += $edges{$e}; } }

    @{$c} = sort { $deg{$b} <=> $deg{$a} } keys %deg;

    foreach my $node (@{$c}) {
        print CC "\t$node"; }
    print CC "\n"; }

close CC;
