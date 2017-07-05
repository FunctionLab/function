#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
use PDF;
use Time::SoFar qw(runtime);
use Graph::Undirected;

my $time;
$time = runtime(); print "\n$time: Reading in networks ...\n";

my($in_net1, $in_net2, $out_file) = @ARGV;
open HH, ">$out_file";

$time = runtime(); print "$time: Indexing exclusion edge list ...\n";

my %ex_edges = (); my (@p, $e);
# open EX, "$in_netex" or die "Can't open $in_netex!";
open EX, "$in_net2" or die "Can't open $in_net2!";
while (<EX>) {
    chomp($_); @p = split '\t', $_;
    $e = join '__', sort($p[0], $p[1]);
    $ex_edges{$e} = $p[2];
}
close EX;

$time = runtime(); print "$time: Indexing subject network ...\n";

my $snet = Graph::Undirected->new;
my %snet_edges = (); my %snet_degree = (); my %snet_neighbors = ();

open SNET, "$in_net1" or die "Can't open $in_net1!";
foreach (<SNET>) {
    chomp($_); @p = split '\t', $_;
    $e = join '__', sort($p[0], $p[1]);
    if((exists $snet_edges{$e}) or ($p[2] <= 0)) { next; }
    if((exists $ex_edges{$e}) and ($ex_edges{$e} > 0)) { next; }

    $snet_edges{$e} = $p[2];
    $snet_degree{$p[0]} += $p[2]; $snet_degree{$p[1]} += $p[2];
    $snet_neighbors{$p[0]}{$p[1]}++; $snet_neighbors{$p[1]}{$p[0]}++;
    $snet->add_weighted_edge($p[0], $p[1], $p[2]);
}
close SNET;

my $tot_ngenes = scalar keys %snet_degree;
my $tot_wgenes = 0;
foreach (keys %snet_degree) {
    $tot_wgenes += (1 - ($snet_degree{$_}/$tot_ngenes)); }

# my $tot_edges = $snet->edges;
# my $edge_prob = 2*$tot_edges/($tot_ngenes*($tot_ngenes-1));

$time = runtime(); print "$time: Indexing weighted subject network ...\n";

my $w1snet = Graph::Undirected->new; # my %w1snet_edge = ();
# my $w2snet = Graph::Undirected->new; # my %w2snet_edge = ();
my %w1snet_degree = ();
my ($w_uv, $p_u, $p_v);
foreach my $e (keys %snet_edges) {
    @p = split '__', $e;
    # $w = $tot_ngenes/sqrt($snet_degree{$p[0]}*$snet_degree{$p[1]});
    # $w = $tot_ngenes/($snet_degree{$p[0]}*$snet_degree{$p[1]});
    $p_u = ($snet_degree{$p[0]}/$tot_ngenes);
    $p_v = ($snet_degree{$p[1]}/$tot_ngenes);
    
    $w_uv = sqrt((1-$p_u)*(1-$p_v));
    $w1snet_degree{$p[0]} += $w_uv; $w1snet_degree{$p[1]} += $w_uv;
    $w1snet->add_weighted_edge($p[0], $p[1], $w_uv);

    # $w_uv = sqrt($p_u*$p_v);
    # $w2snet->add_weighted_edge($p[0], $p[1], $w_uv);
}

$time = runtime(); print "$time: Calculating TOM & shortest-paths ...\n";
print "\t", $tot_ngenes*($tot_ngenes-1)/2, " edges\n";

my ($count, $com_nngb, $com_wngb, $ntom, $wtom, $tempg, @path, $length);
my ($uni_nngb, $uni_wngb, $njac, $wjac, $npval, $wpval);
my ($min_wdeg, $min_ndeg, $min_node, $max_wdeg, $max_ndeg, $max_node);

print HH "Node1\tNode2\tnDeg1\tnDeg2\tProb\twDeg1\twDeg2\tCom.nNgb\tCom.wNgb";
print HH "\tnTOM\twTOM\tnJac\twJac\tnPval\twPval\n";
if($in_net2 eq $in_net1) {
    my @node_array = keys %snet_degree; my ($u, $v);
    for(my $i=0; $i<$#node_array; $i++) {
        $u = $node_array[$i];
        for(my $j=($i+1); $j<=$#node_array; $j++) {
            $v = $node_array[$j];
            $count++; unless($count % 100000) {
                $time = runtime(); print "\t$time: $count ...\n"; }

            $max_wdeg = $w1snet_degree{$u}; $max_node = $u;
            $min_wdeg = $w1snet_degree{$v}; $min_node = $v;
            if($w1snet_degree{$v} > $max_wdeg) {
                $max_wdeg = $w1snet_degree{$v};
                $min_wdeg = $w1snet_degree{$u};
                $max_node = $v;
                $min_node = $u;
            }

            $min_ndeg = $snet_degree{$u}; if($snet_degree{$v} < $min_ndeg) {
                $min_ndeg = $snet_degree{$v}; }
            if($w1snet->has_edge($u, $v)) {
                $max_wdeg -= $w1snet->get_edge_weight($u, $v);
                $min_wdeg -= $w1snet->get_edge_weight($u, $v);
                $min_ndeg--; 
            }

            $com_nngb = 0; $com_wngb = 0;
            foreach my $n (keys %{$snet_neighbors{$u}}) {
                if(exists $snet_neighbors{$v}{$n}) {
                    $com_nngb++;
                    $com_wngb += $w1snet->get_edge_weight($u, $n) * 
                        $w1snet->get_edge_weight($v, $n);
                }
            }

            $ntom = 0; if($min_ndeg > 0) { $ntom = $com_nngb/$min_ndeg; }
            $wtom = 0; if($min_wdeg > 0) { $wtom = $com_wngb/$min_wdeg; }

            if($com_nngb == 0) { next; }
            # $wtom = 0; if($max_wdeg > 0) { $wtom = $com_wngb/$max_wdeg; }
            print HH "$max_node\t$min_node\t$snet_degree{$max_node}\t$snet_degree{$min_node}\t";
            print HH sprintf("%.6g\t%.6g\t%.6g\t%d\t%d\t%.6g\t%.6g\t%.6g\n",
                sqrt(($snet_degree{$max_node}/$tot_ngenes)*($snet_degree{$min_node}/$tot_ngenes)),
                $max_wdeg, $min_wdeg, $min_ndeg, $com_nngb, $com_wngb, $ntom, $wtom);
            # print HH "\t$max_wdeg\t$min_wdeg\t$min_ndeg\t$com_nngb\t$com_wngb\t$ntom\t$wtom\n";
            # print HH "\t", sprintf("%.6g", $max_wdeg), "\t", sprintf("%.6g", $min_wdeg), "\t$min_ndeg";
            # print HH "\t$com_nngb\t", sprintf("%.6g", $com_wngb), "\t";
            # print HH sprintf("%.6g", $ntom), "\t", sprintf("%.6g", $wtom), "\n";
        }
    }
}
else {
    my %qnet_edges = (); my ($u, $v, $w);

    open QNET, "$in_net2" or die "Can't open $in_net2!";
    foreach (<QNET>) {
        chomp($_); ($u, $v, $w) = split '\t', $_;
        $e = join '__', sort($u, $v);

        if((exists $qnet_edges{$e}) or ($u eq $v)) { next; }

        if((exists $snet_degree{$u}) and (exists $snet_degree{$v})) {
            $qnet_edges{$e}++;
            $count++; unless($count % 1000) {
                $time = runtime(); print "\t$time: $count ...\n"; }

            $max_wdeg = $w1snet_degree{$u}; $max_node = $u;
            $min_wdeg = $w1snet_degree{$v}; $min_node = $v;
            if($w1snet_degree{$v} > $max_wdeg) {
                $max_wdeg = $w1snet_degree{$v};
                $min_wdeg = $w1snet_degree{$u};
                $max_node = $v;
                $min_node = $u;
            }

            $min_ndeg = $snet_degree{$u}; $max_ndeg = $snet_degree{$v};
            if($max_ndeg < $min_ndeg) {
                $min_ndeg = $snet_degree{$v}; $max_ndeg = $snet_degree{$u}; }

            if($snet->has_edge($u, $v)) {
                $max_wdeg -= $w1snet->get_edge_weight($u, $v);
                $min_wdeg -= $w1snet->get_edge_weight($u, $v);
                $max_ndeg--; $min_ndeg--; 
            }

            $com_nngb = $com_wngb = $uni_nngb = $uni_wngb = 0;
            foreach my $n (keys %{$snet_neighbors{$u}}) {
                $uni_nngb++;
                if(exists $snet_neighbors{$v}{$n}) {
                    $com_nngb++;
                    $com_wngb += $w1snet->get_edge_weight($u, $n) * 
                        $w1snet->get_edge_weight($v, $n);

                    $uni_wngb += (sort ($w1snet->get_edge_weight($u, $n),
                        $w1snet->get_edge_weight($v, $n)))[1];
                }
                else {
                    $uni_wngb += $w1snet->get_edge_weight($u, $n); }
            }

            foreach my $n (keys %{$snet_neighbors{$v}}) {
                if(exists $snet_neighbors{$u}{$n}) { next; }

                $uni_nngb++;
                $uni_wngb += $w1snet->get_edge_weight($v, $n);
            }

            $ntom = 0; if($min_ndeg > 0) { $ntom = $com_nngb/$min_ndeg; }
            $wtom = 0; if($min_wdeg > 0) { $wtom = $com_wngb/$min_wdeg; }

            $njac = $com_nngb/$uni_nngb;
            $wjac = $com_wngb/$uni_wngb;

            $npval = 1; if ($com_nngb > 0) {
                $npval = hypergeometric_tail($tot_ngenes, $min_ndeg,
                    $max_ndeg, $com_nngb); }
            $npval = (-1)*log($npval)/log(10);

            # $min_wdeg = int($min_wdeg + 0.5);
            # $max_wdeg = int($max_wdeg + 0.5);

            $wpval = 1; if(int($com_wngb + 0.5) > 0) {
                $wpval = hypergeometric_tail(int($tot_wgenes + 0.5), int($min_wdeg + 0.5),
                    int($max_wdeg + 0.5), int($com_wngb + 0.5)); }
            $wpval = (-1)*log($wpval)/log(10);

            # if($com_nngb == 0) { next; }

            # if(exists $snet_edges{$e}) {
                # $tempg = $w2snet->copy_graph;
                #     $tempg = $snet->copy_graph;
                # $tempg->delete_edge($u, $v);
                # @path = $tempg->SP_Dijkstra($u, $v);
                # }
            # else { @path = $w2snet->SP_Dijkstra($u, $v); }
            # else { @path = $snet->SP_Dijkstra($u, $v); }
            # $length = $#path;
            # $length = 0;

            print HH "$max_node\t$min_node\t$snet_degree{$max_node}\t$snet_degree{$min_node}\t";
            print HH sprintf("%.6g\t%.6g\t%.6g\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g",
                sqrt((1 - $snet_degree{$max_node}/$tot_ngenes)*(1 - $snet_degree{$min_node}/$tot_ngenes)),
                $max_wdeg, $min_wdeg, $com_nngb, $com_wngb, $ntom, $wtom, $njac, $wjac, $npval, $wpval), "\n";
            # "\t$length\n";
        }
    }
    close QNET;
}

$time = runtime(); print "$time: DONE\n\n";

close HH;
