#!/usr/bin/perl

# ================================================
# Name : pick-uniq.pl
# Purpose : Pick 'uniq' nodes from an overlap graph
# Created : 10-03-2015
# Last Modified : Wed 16 Sep 2015 05:51:28 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Graph::Undirected;


my $iterm = $ARGV[0];   # terms and enrichemnt scores
my $iovlp = $ARGV[1];   # highly overlapping term pairs based on ovlp; output from gsea_hypergeometric: term1 term2 ovlp
my $isubgs = $ARGV[2];  # list of initially tested terms
my $irel = $ARGV[3];    # ontological relationships from propagate_annotations
my $otab = $ARGV[4];    # output

#my $icut_ovp = 0.9;
#my $icut_jac = 0.5;
my $icut_ovp = $ARGV[5];    # Cutoff for overlap
my $icut_jac = $ARGV[6];    # Cutoff for jaccard


my %sub_gs = ();
open SGS, "$isubgs" or die "Can't open $isubgs!";
while (<SGS>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $sub_gs{$p[0]}++; }
close SGS;


my $ovp_graph = Graph::Undirected->new;
my %node_score = my %node_desc = my %node_size = ();
open GH, "$iterm" or die "Can't open $iterm!";
while (<GH>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    unless(exists $sub_gs{$p[0]}) { next; }
    $ovp_graph->add_vertex($p[0]);
    $node_score{$p[0]} = $p[$#p];
    $node_desc{$p[0]} = $p[2];
    $node_size{$p[0]} = $p[1]; }
close GH;


my $rel_graph = Graph::Undirected->new;
open REL, "$irel" or die "Can't open $irel!";
while (<REL>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    unless(exists $sub_gs{$p[0]}) { next; }

    my $ch = shift @p; shift @p;
    foreach my $pa (@p) {
        unless(exists $sub_gs{$pa}) { next; }
        if($pa eq $ch) { next; }
        $rel_graph->add_edge($ch, $pa); } }
close REL;

my @rel_cc = $rel_graph->connected_components();

my %rel_pairs = ();
foreach my $rc (@rel_cc) {
    for(my $i=0; $i<$#{$rc}; $i++) {
        for(my $j=($i+1); $j<=$#{$rc}; $j++) {
            my $e = join '__', sort(${$rc}[$i], ${$rc}[$j]);
            $rel_pairs{$e}++; } } }


my %node_wdeg = my %node_nei = my %jac_pairs = ();

open OVP, "$iovlp" or die "Can't open $iovlp!";
while (<OVP>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $gs1 = $p[0];
    my $gs2 = $p[3];
    my $ovp = $p[7];
    my $jac = $p[8];

    if($ovp < $icut_ovp) { next; }
    unless((exists $node_score{$gs1}) and (exists $node_score{$gs2})) {
        next; }
    
    $ovp_graph->add_edge($gs1, $gs2);

    $node_wdeg{$gs1} += $ovp;
    $node_wdeg{$gs2} += $ovp;

    $node_nei{$gs1}{$gs2} = $ovp;
    $node_nei{$gs2}{$gs1} = $ovp;

    if($jac < $icut_jac) { next; }
    my $e = join '__', sort($gs1, $gs2);

    $jac_pairs{$e}++; }
close OVP;


my @cc = $ovp_graph->connected_components();


my $ncc = 0; my %cc_maxs = my %cc_nodes = ();
foreach my $c (@cc) {
    $ncc++; my $ovpcc_tag = sprintf("%04d", $ncc);
    my $n = scalar @$c;

    if($n == 1) {
        my $top_node = ${$c}[0];
        push(@{$cc_nodes{$ovpcc_tag}}, $top_node);
        $cc_maxs{$ovpcc_tag} = $node_score{$top_node}; }

    else {
        my %rem_nodes = ();
        map { $rem_nodes{$_}++; } @$c;
        my $nrem; my $itr = 0;

        do {
            $itr++;
            my @node_srt = sort { $node_wdeg{$b} <=> $node_wdeg{$a} ||
                $node_score{$b} <=> $node_score{$a} } keys %rem_nodes;

            my $top_node = $node_srt[0];
            push(@{$cc_nodes{$ovpcc_tag}}, $top_node);
            if(exists $cc_maxs{$ovpcc_tag}) {
                if($node_score{$top_node} > $cc_maxs{$ovpcc_tag}) {
                    $cc_maxs{$ovpcc_tag} = $node_score{$top_node}; } }
            else { $cc_maxs{$ovpcc_tag} = $node_score{$top_node}; }
            delete $rem_nodes{$top_node};

            foreach my $nei (keys %{$node_nei{$top_node}}) {
                $node_wdeg{$nei} -= $node_nei{$top_node}{$nei};
                delete $node_nei{$nei}{$top_node}; }
            delete $node_nei{$top_node};

            $nrem = scalar keys %rem_nodes;
        } until ($nrem == 0); } }
print "\n$ncc connected components.\n\n";


my %ontcc_maxs = my %jaccc_maxs = ();
open TAB, ">$otab";
foreach my $ovpcc_tag (sort {$cc_maxs{$b} <=> $cc_maxs{$a}} keys %cc_maxs) {
    my %ont_jac_nodes = ();
    my $n = scalar @{$cc_nodes{$ovpcc_tag}};
    print "$ovpcc_tag\t$n\n";

    my $ccrel_graph = Graph::Undirected->new;
    
    my %cc_node_rank = ();
    for(my $i=0; $i<=$#{$cc_nodes{$ovpcc_tag}}; $i++) {
        $ccrel_graph->add_vertex(${$cc_nodes{$ovpcc_tag}}[$i]);
        $cc_node_rank{${$cc_nodes{$ovpcc_tag}}[$i]} = $i; }

    for(my $i=0; $i<$#{$cc_nodes{$ovpcc_tag}}; $i++) {
        my $node1 = ${$cc_nodes{$ovpcc_tag}}[$i];

        for(my $j=($i+1); $j<=$#{$cc_nodes{$ovpcc_tag}}; $j++) {
            my $node2 = ${$cc_nodes{$ovpcc_tag}}[$j];

            my $e = join '__', sort($node1, $node2);
            if(exists $rel_pairs{$e}) {
                $ccrel_graph->add_edge($node1, $node2); } } }

    my @ont_cc = $ccrel_graph->connected_components();

    my $nont_cc = 0;
    foreach my $sc (@ont_cc) {
        $nont_cc++; my $ontcc_tag = $ovpcc_tag.'.'.sprintf("%04d", $nont_cc);

        #my @anodes = sort { $cc_node_rank{$a} <=> $cc_node_rank{$b}} @$sc;
        #print "\t$ontcc_tag\t", scalar @anodes, "\n";

        my $ccjac_graph = Graph::Undirected->new;
        #foreach my $node (@anodes) {
        foreach my $node (@$sc) {
            if(exists $ontcc_maxs{$ontcc_tag}) {
                if($node_score{$node} > $ontcc_maxs{$ontcc_tag}) {
                    $ontcc_maxs{$ontcc_tag} = $node_score{$node}; } }
            else {
                $ontcc_maxs{$ontcc_tag} = $node_score{$node}; }
            $ccjac_graph->add_vertex($node); }

        for(my $i=0; $i<$#{$sc}; $i++) {
            for(my $j=($i+1); $j<=$#{$sc}; $j++) {
                my $e = join '__', sort(${$sc}[$i], ${$sc}[$j]);
                if(exists $jac_pairs{$e}) {
                    $ccjac_graph->add_edge(${$sc}[$i], ${$sc}[$j]); } } }
                    #my $n1 = $anodes[$i]; my $n2 = $anodes[$j];
                    #print ">>$ovpcc_tag\t$ontcc_tag\t$n1|$n2\t";
                    #print "$node_desc{$n1}|$node_desc{$n2}\t";
                    #print "$node_size{$n1}|$node_size{$n2}\t";
                    #print "$node_score{$n1}|$node_score{$n2}\n"; } } }

        my @jac_cc = $ccjac_graph->connected_components();

        my $njac_cc = 0;
        foreach my $jc (@jac_cc) {
            $njac_cc++; my $jaccc_tag = $ontcc_tag.'.'.sprintf("%04d", $njac_cc);
            #my @anodes2 = sort { $cc_node_rank{$a} <=> $cc_node_rank{$b}} @$jc;
            #print "\t\t$jaccc_tag\t", scalar @anodes2, "\n";

            #foreach my $node (@anodes2) {
            foreach my $node (@$jc) {
                if(exists $jaccc_maxs{$jaccc_tag}) {
                    if($node_score{$node} > $jaccc_maxs{$jaccc_tag}) {
                        $jaccc_maxs{$jaccc_tag} = $node_score{$node}; } }
                else {
                    $jaccc_maxs{$jaccc_tag} = $node_score{$node}; }
                $ont_jac_nodes{$ontcc_tag}{$jaccc_tag}{$node}++; } }
                #print TAB "$ovpcc_tag|$ontcc_tag|$jaccc_tag\t$node\t";
                #print TAB "$node_desc{$node}\t$node_size{$node}\t$node_score{$node}\n"; }
            #print TAB "----------------------------------------------------------------------\n"; }

            #print TAB "======================================================================\n";
    }

    print TAB "$ovpcc_tag ($n)\n";

    my @aontcc = sort { $ontcc_maxs{$b} <=> $ontcc_maxs{$a} } keys %ont_jac_nodes;
    foreach my $ontcc_tag (@aontcc) {
        (my $tago = $ontcc_tag) =~ s/^$ovpcc_tag\.//g;
        print TAB "- $tago ------------------------------------------------------------------------\n";

        my @ajaccc = sort { $jaccc_maxs{$b} <=> $jaccc_maxs{$a} } keys %{$ont_jac_nodes{$ontcc_tag}};
        foreach my $jaccc_tag (@ajaccc) {
            (my $tagj = $jaccc_tag) =~ s/^$ontcc_tag\.//g;
            #print TAB "\t\t$jaccc_tag\n";

            my @anodes = sort { $cc_node_rank{$a} <=> $cc_node_rank{$b}} keys %{$ont_jac_nodes{$ontcc_tag}{$jaccc_tag}} ;
            my $njcc = scalar @anodes;
            foreach my $node (@anodes) {
                #print TAB "-- $tagj (", sprintf("%02d", $njcc),")";
                my $dots = '.' x $njcc; $njcc--;
                print TAB "-- $tagj $dots ";
                print TAB "$node\t$node_desc{$node}\t$node_size{$node}\t$node_score{$node}\n";
            }
        }
    }

    print TAB "================================================================================\n"; }
close TAB;

