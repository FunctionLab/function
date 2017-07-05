#!/usr/bin/perl
use Time::SoFar qw( runtime runinterval figuretimes );

$f1=$ARGV[0];		# Net1: SIF or DAT format
$f2=$ARGV[1];		# Net2: SIF or DAT format
$anal=$ARGV[2];		# Scoring: 'min', 'max', 'sum', 'diff', 'prod', 'avg', 'gmean'
$norm=$ARGV[3];		# Normalization: '0to1', 'mean-center', 'colz', 'none'
$pop=$ARGV[4];		# Print option - Union or Intersection: 'u' or 'i'
chomp($f3=$ARGV[5]);	# Outfile: Combined network

open FH, "$f1" or die "Can't open $f1!"; chomp(my @f=<FH>); close FH;
open GH, "$f2" or die "Can't open $f2!"; chomp(my @g=<GH>); close GH;
open HH, ">$f3";

$elapsed = runtime(); print "\n$elapsed: Reading and indexing both networks ... \n";

%net1genes=(); %net1edges=(); $maxnet1=0;
foreach(@f) {
	$_=~s///g; @p=split("\t",$_);

    $edge = join '__', sort($p[0], $p[1]);

	$net1genes{$p[0]}++; $net1genes{$p[1]}++;
	$net1edges{$edge}=$p[2];

	if($p[2]>$maxnet1) { $maxnet1=$p[2]; }
}

%net2genes=(); %net2edges=(); $maxnet2=0;
foreach(@g) {
	$_=~s///g; @p=split("\t",$_);

    $edge = join '__', sort($p[0], $p[1]);

	$net2genes{$p[0]}++; $net2genes{$p[1]}++;
	$net2edges{$edge}=$p[2];

	if($p[2]>$maxnet2) { $maxnet2=$p[2]; }
}

$numgenesnet1=0; $numedgesnet1=0;
$numgenesnet2=0; $numedgesnet2=0;
%commongenes=(); %uniongenes=();

foreach $gene (keys %net1genes) {
	$numgenesnet1++;
	$uniongenes{$gene}++;

	if(exists $net2genes{$gene}) { $commongenes{$gene}++; }
}

foreach $gene (keys %net2genes) {
	$numgenesnet2++;
	$uniongenes{$gene}++;
}

$elapsed = runtime(); print "\n$elapsed: Comparing the networks ... \n";

%commonedges=(); %unionedges=();
foreach $edge (keys %net1edges) {
	$numedgesnet1++; $net1score=$net1edges{$edge};
	if($norm eq '0to1') { $net1score=($net1score/$maxnet1); }

	if(exists $net2edges{$edge}) {
		$commonedges{$edge}++; $net2score=$net2edges{$edge};
		if($norm eq '0to1') { $net2score=($net2score/$maxnet2); }

		$maxscore=0; $minscore=0;
		if($net1score<$net2score) {
			$minscore=$net1score;
			$maxscore=$net2score;
		}
		else {
			$minscore=$net2score;
			$maxscore=$net1score;
		}

		$sign1=0; if($net1score >= 0) { $sign1=1; } else { $sign1=-1; }
		$sign2=0; if($net2score >= 0) { $sign2=1; } else { $sign2=-1; }

		if($anal eq 'min') { $unionedges{$edge}=$minscore; }
		if($anal eq 'max') { $unionedges{$edge}=$maxscore; }
		if($anal eq 'avg') { $avgscore=(($net1score+$net2score)/2); $unionedges{$edge}=$avgscore; }
		if($anal eq 'gmean') { $gmeanscore=($sign1*$sign2)*(sqrt(abs($net1score)*abs($net2score))); $unionedges{$edge}=$gmeanscore; }
		if($anal eq 'sum') { $sumscore=($net1score+$net2score); $unionedges{$edge}=$sumscore; }
		if($anal eq 'prod') { $prodscore=($net1score*$net2score); $unionedges{$edge}=$prodscore; }
		if($anal eq 'diff') { $diffscore=($net1score-$net2score); $unionedges{$edge}=$diffscore; }

	}
	else {
		$unionedges{$edge}=$net1score;
	}
}

foreach $edge (keys %net2edges) {
	$numedgesnet2++;

	unless(exists $net1edges{$edge}) {
		$unionedges{$edge}=$net2edges{$edge};
	}
}

$numcommongenes=0; $numcommonedges=0;
$numuniongenes=0; $numunionedges=0;

foreach (keys %commongenes) { $numcommongenes++; }
foreach (keys %uniongenes) { $numuniongenes++; }

$elapsed = runtime(); print "\n$elapsed: Printing edges of choice ... \n";

%printgenes=();
if($pop eq 'u') {
	foreach $edge (sort {$unionedges{$b} <=> $unionedges{$a}} keys %unionedges) {
		$numunionedges++;
		@p=split("__",$edge);

		print HH "$p[0]\t$p[1]\t$unionedges{$edge}\n";

		$printgenes{$p[0]}++; $printgenes{$p[1]}++;
	}
}
elsif($pop eq 'i') {
	foreach $edge (sort {$commonedges{$b} <=> $commonedges{$a}} keys %commonedges) {
		$numcommonedges++;
		@p=split("__",$edge);

		$n1=sprintf("%.6g", $net1edges{$edge});
        $n2=sprintf("%.6g", $net2edges{$edge});
        $u=sprintf("%.6g", $unionedges{$edge});
		
        print HH "$p[0]\t$p[1]\t$n1\t$n2\t$u\n";
	
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

$elapsed = runtime(); print "$elapsed: DONE \n\n";

close HH;

