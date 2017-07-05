#!/usr/bin/perl
use strict;
use warnings;

my ($ilist, $omat) = @ARGV;

my @col = qw(lin exp1 exp2 exp3);
my %pwr = qw(lin 0 exp1 4 exp2 3.5 exp3 3);


open MAT, ">$omat";
print MAT "gene";
foreach my $c (@col) {
    print MAT "\t$c"; }
print MAT "\n";


my %gene_score = ();
open LI, "$ilist" or die "Can't open $ilist!";
while (<LI>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $gene_score{$p[0]} = $p[1]; }
close LI;


my $N = scalar keys %gene_score;
my @sortg = sort {$gene_score{$b} <=> $gene_score{$a}} keys %gene_score;

my $r = 0;
foreach my $g (@sortg) {
    print MAT "$g"; $r++;
    foreach my $c (@col) {
        # print "$c\t$pwr{$c}\t",(10**$pwr{$c}),"\t",(1/20**$pwr{$c}),"\n";
        if($pwr{$c} == 0) {
            printf MAT "\t%.6g", ( ( $N - $r + 1 ) / $N ); }
        else {
            printf MAT "\t%.6g", (exp( -2*($r-1)/10**$pwr{$c} )) ; } }
    # exit;
    print MAT "\n"; }
close MAT;
