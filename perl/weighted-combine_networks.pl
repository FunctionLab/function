#!/usr/bin/perl
use Time::SoFar qw( runtime runinterval figuretimes );

$f1=$ARGV[0];		# Net1: DAT format; Edge weights in condition 1
$f2=$ARGV[1];		# Net2: DAT format; Edge weights in condition 2
$f3=$ARGV[2];		# Net3: DAT format; This file contains the original edge weights
$mode=$ARGV[3];		# Scoring: 'min', 'max', 'sum', 'diff', 'prod', 'avg', 'gmean', 'weighted-sum', 'weighted-diff'
$norm=$ARGV[4];		# Normalization: '0to1', 'mean-center', 'colz', 'none'
$pop=$ARGV[5];		# Print option - Union or Intersection: 'u' or 'i'
chomp($f4=$ARGV[6]);	# Outfile: Combined network

open N1H,"$f1";
chomp(@net1=<N1H>);
close N1H;
open N2H,"$f2";
chomp(@net2=<N2H>);
close N2H;
open N3H,"$f3";
chomp(@net3=<N3H>);
close N3H;

open HH,">$f4";

my (@p, $e, $time);
$time = runtime(); print "\n$time: Reading and indexing both networks ... \n";

my %net1g = my %net1e = (); my $net1_maxe=0;
open NET, "$inet1" or die "Can't open $inet1!";
while (<NET>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $net1g{$p[0]}++; $net1g{$p[1]}++;
    $e = join '__', sort($p[0], $p[1]);
    $net1e{$e} = $p[2];

    if($p[2] > $net1_maxe) { $net1_maxe = $p[2]; } }
close NET;

my %net2g = my %net2e = (); my $net2_maxe=0;
open NET, "$inet2" or die "Can't open $inet2!";
while (<NET>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $net2g{$p[0]}++; $net2g{$p[1]}++;
    $e = join '__', sort($p[0], $p[1]);
    $net2e{$e} = $p[2];

    if($p[2] > $net2_maxe) { $net2_maxe = $p[2]; } }
close NET;

my %net3g = my %net3e = (); my $net3_maxe=0;
open NET, "$inet3" or die "Can't open $inet3!";
while (<NET>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $net3g{$p[0]}++; $net3g{$p[1]}++;
    $e = join '__', sort($p[0], $p[1]);
    $net3e{$e} = $p[2];

    if($p[2] > $net3_maxe) { $net3_maxe = $p[2]; } }
close NET;

$numgenesnet1=0; $numedgesnet1=0;
$numgenesnet2=0; $numedgesnet2=0;
$numgenesnet3=0; $numedgesnet3=0;

%commongenes=(); %uniongenes=();

foreach $gene (keys %net1g)
{
	$numgenesnet1++;
	$uniongenes{$gene}++;

	if(exists $net2g{$gene})
	{
		$commongenes{$gene}++;
	}
}

foreach $gene (keys %net2g)
{
	$numgenesnet2++;
	$uniongenes{$gene}++;
}

$time = runtime(); print "\n$time: Comparing the networks ... \n";

%commonedges=(); %unionedges=();
foreach $edge (keys %net1e)
{
	$numedgesnet1++; $net1score=$net1e{$edge};
	if($norm eq '0to1') { $net1score=($net1score/$net1_maxe); }

	if((exists $net2e{$edge})&&(exists $net3edges{$edge}))
	{
		$commonedges{$edge}++; $net2score=$net2e{$edge}; $weight=$net3edges{$edge};
		if($norm eq '0to1') { $net2score=($net2score/$net2_maxe); }

		$maxscore=0; $minscore=0;
		if($net1score<$net2score)
		{
			$minscore=$net1score;
			$maxscore=$net2score;
		}
		else
		{
			$minscore=$net2score;
			$maxscore=$net1score;
		}

		$sign1=0; if($net1score >= 0) { $sign1=1; } else { $sign1=-1; }
		$sign2=0; if($net2score >= 0) { $sign2=1; } else { $sign2=-1; }

		if($mode eq 'min') { $unionedges{$edge}=$minscore; }
		if($mode eq 'max') { $unionedges{$edge}=$maxscore; }
		if($mode eq 'avg') { $avgscore=(($net1score+$net2score)/2); $unionedges{$edge}=$avgscore; }
		if($mode eq 'gmean') { $gmeanscore=($sign1*$sign2)*(sqrt(abs($net1score)*abs($net2score))); $unionedges{$edge}=$gmeanscore; }
		if($mode eq 'sum') { $sumscore=($net1score+$net2score); $unionedges{$edge}=$sumscore; }
		if($mode eq 'prod') { $prodscore=($net1score*$net2score); $unionedges{$edge}=$prodscore; }
		if($mode eq 'diff') { $diffscore=($net1score-$net2score); $unionedges{$edge}=$diffscore; }
		if($mode eq 'weighted-sum') { $wsumscore=$weight*($net1score+$net2score); $unionedges{$edge}=$wsumscore; }
		if($mode eq 'weighted-diff') { $wdiffscore=$weight*($net1score-$net2score); $unionedges{$edge}=$wdiffscore; }
	}
	else
	{
		$unionedges{$edge}=$net1score;
	}
}

foreach $edge (keys %net2e)
{
	$numedgesnet2++;

	unless(exists $net1e{$edge})
	{
		$unionedges{$edge}=$net2e{$edge};
	}
}

$numcommongenes=0; $numcommonedges=0;
$numuniongenes=0; $numunionedges=0;

foreach (keys %commongenes) { $numcommongenes++; }
foreach (keys %uniongenes) { $numuniongenes++; }

$time = runtime(); print "\n$time: Printing edges of choice ... \n";

%printgenes=();
if($pop eq 'u')
{
	foreach $edge (sort {$unionedges{$b} <=> $unionedges{$a}} keys %unionedges)
	{
		$numunionedges++;

		@p=();
		@p=split("__",$edge);

		print HH "$p[0]\t$p[1]\t$unionedges{$edge}\n";

		$printgenes{$p[0]}++; $printgenes{$p[1]}++;
	}
}

elsif($pop eq 'i')
{
	foreach $edge (sort {$commonedges{$b} <=> $commonedges{$a}} keys %commonedges)
	{
		$numcommonedges++;

		@p=();
		@p=split("__",$edge);

		$w=sprintf("%.6g", $net3edges{$edge});
		$n1=sprintf("%.6g", $net1e{$edge}); $n2=sprintf("%.6g", $net2e{$edge});
		$u=sprintf("%.6g", $unionedges{$edge});
		print HH "$p[0]\t$p[1]\t$w\t$n1\t$n2\t$u\n";
	
		$printgenes{$p[0]}++; $printgenes{$p[1]}++;
	}
}

$numprintgenes=0; $numprintgenes=scalar(keys %printgenes);

print "\nNet1: $f1\nNo. of genes: $numgenesnet1\nNo. of edges: $numedgesnet1\n";
print "\nNet2: $f2\nNo. of genes: $numgenesnet2\nNo. of edges: $numedgesnet2\n";

print "\nCommon\nNo. of genes: $numcommongenes\nNo. of edges: $numcommonedges\n";
print "\nUnion\nNo. of genes: $numuniongenes\nNo. of edges: $numunionedges\n\n";

if($pop eq 'u') { $printtag='Union'; } elsif($pop eq 'i') { $printtag='Intersection'; }
print "Printing $printtag containing $numprintgenes\n\n";

$time = runtime(); print "$time: DONE \n\n";

close HH;

