#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my($f1, $f2, $multmap, $scorecalc, $f3) = @ARGV;
my $directed = 1;

# $f1: Network in DAT format
# $f2: ID mapping: <id1> <id2>
# $multmap: Allow one-to-many mapping? 'yes' or 'no'
# $scorecalc: Choice of how to find the score of an edge given multiple values when mapped: 'max', 'min', or 'mean'
# $f3: Outfile; Network in new ID space

open FH, "$f1" or die "Can't open $f1!"; chomp(my @net=<FH>); close FH;
open GH, "$f2" or die "Can't open $f2!"; chomp(my @map=<GH>); close GH;
open HH, ">$f3";

my %idmap=(); my @p;
foreach my $id (@map) {
	$id=~s///g; @p=split("\t",$id);

	for(my $j=1; $j<=$#p; $j++) {
		push(@{$idmap{uc($p[0])}}, uc($p[$j]));
	}
}

my %new_edges = (); my ($tot_oldedges, $par) = (0, 0); my $edge;
foreach my $oe (@net) {
	$tot_oldedges++;
	$oe=~s///g; $oe = uc($oe); @p=split("\t",$oe);

	if((exists $idmap{$p[0]}) and (exists $idmap{$p[1]})) {
		$par=0;
		if($multmap eq 'no') {
			if(($#{$idmap{$p[0]}} eq 0) and ($#{$idmap{$p[1]}} eq 0)) {
				$par=1;
			}
		}
		else { $par=1; }

		if($par eq 1) {
            foreach my $from (@{$idmap{$p[0]}}) {
                foreach my $to (@{$idmap{$p[1]}}) {
                    if($directed == 0) { $edge = join '__', sort($from, $to); }
                    else { $edge = join '__', ($from, $to); }
					push(@{$new_edges{$edge}}, $p[2]);
                }
            }
		}
	}
}

my ($tot_newedges, $num_newedges_wmultscores) = (0, 0);
my @sorted_edge_scores = (); my ($sum_edge_score, $edge_score) = (0, 0);
foreach my $ne (keys %new_edges) {
	$tot_newedges++;

	if($#{$new_edges{$ne}} > 0) {
		$num_newedges_wmultscores++;
        @sorted_edge_scores = sort {$a <=> $b} @{$new_edges{$ne}};

		if($scorecalc eq 'max') {
            $edge_score = $sorted_edge_scores[$#sorted_edge_scores]; }
		elsif($scorecalc eq 'min') {
            $edge_score = $sorted_edge_scores[0]; }
		else {
			$sum_edge_score = 0;
			foreach my $s (@{$new_edges{$ne}}) {
				$sum_edge_score += $s;
			}

			$edge_score = ($sum_edge_score/(scalar(@{$new_edges{$ne}})));
		}
	}
	else { $edge_score=${$new_edges{$ne}}[0]; }

	@p = split("__",$ne);
	print HH "$p[0]\t$p[1]\t$edge_score\n";
}

print "\nTot. old edges: $tot_oldedges\n";
print "Tot. new edges: $tot_newedges\n";
print "No. new edges with multiple scores: $num_newedges_wmultscores\n\n";

close HH;

