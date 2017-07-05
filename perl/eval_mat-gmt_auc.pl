#!/usr/bin/perl
use strict;
use warnings;

sub dcheck;

my ($imat, $igmt, $oeval) = @ARGV;


my %gene_truel = my %all_gmt_genes = ();
open GMT, "$igmt" or die "Can't open $igmt!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $gs = shift @p; shift @p;
    foreach my $g (@p) {
        $all_gmt_genes{$g}++;
        $gene_truel{$gs}{$g} = 1; } }
close GMT;


open MH, "$imat" or die "Can't open $imat!";
my @tissues = my %gene_preds = my %all_pred_genes = ();
while (<MH>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    if($_ =~ /^Gene/) {
        @tissues = @p; next; }

    unless(exists $all_gmt_genes{$p[0]}) { next; }
    $all_pred_genes{$p[0]}++;

    for(my $i=1; $i<=$#p; $i++) {
        if($p[$i] =~ /inf/) { next; }
        $gene_preds{$tissues[$i]}{$p[0]} = abs($p[$i]); } }
close MH;


open EV, ">$oeval"; print EV "Tissue";

my @gs_array = ();
foreach my $gs (sort keys %gene_truel) {
    foreach my $g (keys %all_gmt_genes) {
        if(exists $gene_truel{$gs}{$g}) { next; }
        unless(exists $all_pred_genes{$g}) { next; }
        $gene_truel{$gs}{$g} = -1; }

    foreach my $g (keys %{$gene_truel{$gs}}) {
        if(exists $all_pred_genes{$g}) { next; }
        delete $gene_truel{$gs}{$g}; }

    my $P = grep { $_ == 1 } values %{$gene_truel{$gs}};
    my $N = grep { $_ == -1 } values %{$gene_truel{$gs}};
    if(($P == 0) or ($N == 0)) { delete $gene_truel{$gs}; next; }

    push(@gs_array, $gs);
    print "$gs\t$P\t$N\n";
    print EV "\t$gs"; }
print EV "\n";


print "\n"; shift @tissues;
foreach my $t (@tissues) {
    print EV "$t"; print "$t\n";
    
    foreach my $gs (@gs_array) {
        my $auc = dcheck($gene_truel{$gs}, $gene_preds{$t});
        printf EV "\t%.5g", $auc; }

    print EV "\n"; }
print "\n";


# Subroutines
sub dcheck {
    my $truel = shift;
    my $preds = shift;

    my $P = my $N = 0;
    foreach (keys %$truel) {
        unless(exists $preds->{$_}) { next; }
        if($truel->{$_} == 1) { $P++; }
        elsif($truel->{$_} == -1) { $N++; } }

    if(($P == 0) || ($N == 0)) { return; }

    my $tp = my $fp = my $tn = my $fn = my $auc = 0;

    foreach (sort {$preds->{$b} <=> $preds->{$a}} keys %{$preds}) {
        if($truel->{$_} == 1) {
            $tp++; }

        elsif($truel->{$_} == -1) {
            $fp++;
            $auc += $tp; }

        else { next; }

        $tn = $N - $fp; $fn = $P - $tp; }

    if(($tp == 0) or ($fp == 0)) {
        print STDERR "warning: Too few +ve true labels or -ve true labels\n";
        $auc = 0; }
    else {
        $auc = $auc / $tp / $fp; }

    return $auc; }

sub rank {
    my $href = shift;
    my @asort = sort {$a <=> $b} values %$href;
    
    my $i = 1; my %hsort = my %rep = ();
    foreach (@asort) {
        $rep{$_}++;
        if(exists $hsort{$_}) {
            $hsort{$_} = $hsort{$_} + (($i - $hsort{$_}) / $rep{$_}); }
        else { $hsort{$_} = $i; }
        $i++; }

    my %rank = ();
    foreach (keys %$href) { $rank{$_} = $hsort{$href->{$_}}; }
    return \%rank; }

sub mean_sd {
    my $href = shift;

    my $om = my $m = my $s = my $k = 0;
    foreach (keys %$href) {
        $k++;
        my $v = $href->{$_};
        $m = $om + ($v - $om) / $k;
        $s = $s + ($v - $om)*($v - $m);
        $om = $m; }

    $s = sqrt($s / ($k-1));

    return ($m, $s); }

sub cor {
    my $hrefx = shift;
    my $hrefy = shift;

    my $rankx = rank($hrefx);
    my $ranky = rank($hrefy);

    my ($mx, $sx) = mean_sd($rankx);
    my ($my, $sy) = mean_sd($ranky);

    my $n = my $cor = 0;
    foreach my $g (keys %$rankx) {
        $n++;
        $cor += (($rankx->{$g} - $mx) / $sx) * (($ranky->{$g} - $my) / $sy); }
    $cor /= ($n-1);

    return $cor; }

sub rbo {
    my $hrefx = shift;
    my $hrefy = shift;
    my $p = shift;

    my @rankx = sort {$hrefx->{$b} <=> $hrefx->{$a}} keys %$hrefx;
    my @ranky = sort {$hrefy->{$b} <=> $hrefy->{$a}} keys %$hrefy;

    my %rx = my %ry = ();
    my $rbo = 0; # my $p = 0.9;

    for(my $i=0; $i<=$#rankx; $i++) {
        $rx{$rankx[$i]}++;
        $ry{$ranky[$i]}++;

        my $c = grep { exists $rx{$_} } keys %ry;
        
        $rbo += ($c / ($i+1)) * ($p ** $i); }

    $rbo *= (1 - $p);
    return $rbo; }


