#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(shuffle);

my $in_ssv = $ARGV[0];
open SSV, "$in_ssv" or die "Can't open $in_ssv!";

my %dat = (); my @lab = (); my $i = 0; my @p;
while(<SSV>) {
    @p = split ' ', $_;

    push(@lab, shift(@p));
    push(@{$dat{$i}}, @p};

    $i++;
}
close SSV;

(my $out_ssv = $in_ssv) =~ s/.ssv$/.rand-lab.ssv/g;
if($in_ssv != /.ssv$/) { $out_ssv .= '.rand-lab'; }
open RND, ">$out_ssv";

my @rnd_idx = shuffle(keys %dat);
$i = 0; do {
    print RND "$lab[$i] @{$dat{$rnd_idx[$i]}}\n";
    $i++;
} while ($i <= $#rnd_idx);

close RND;
