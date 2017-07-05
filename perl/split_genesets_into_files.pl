#!/usr/bin/perl

$f1=''; $f2=''; $f3='';
chomp($f1=$ARGV[0]);	# GeneSets: <Geneset1> <Gene>...
$f2=$ARGV[1];		# Genes to include
chomp($f3=$ARGV[2]);	# GeneSets to include: <GeneSet> <Desc>

open FH,"$f1"; chomp(@f=<FH>); close FH;

unless($f2 eq '') {
	open GH,"$f2"; chomp(@g=<GH>); close GH;

	%genes=();
	foreach(@g) {
		$_=~s///g;
		@p=(); @p=split("\t",$_);
		$genes{$p[0]}++;
	}
}

%gs_numgenes = ();
unless($f3 eq '') {
	open JH,"$f3"; chomp(@s=<JH>); close JH;

	%geneset=();
	foreach(@s) {
		$_=~s///g;
		@p=(); @p=split("\t",$_);

		$gs_numgenes{$p[0]}=0;
		$gsname=$p[1]; $gsname=~s/ /_/g;
		$geneset{$p[0]}=$gsname;
	}
}

%gs_genes=(); %gsgene_idx=();
foreach(@f) {
	$_=~s///g;
	@p=(); @p=split("\t",$_);

	unless(exists $gsgene_idx{$_}) {
		$gsgene_idx{$_}++;
		if($f2 ne '') {
			if(exists $genes{$p[1]}) {
				if($f3 ne '') {
					if(exists $geneset{$p[0]}) {
						push(@{$gs_genes{$p[0]}}, $p[1]);
						$gs_numgenes{$p[0]}++;
					}
				}
				else {
					push(@{$gs_genes{$p[0]}}, $p[1]);
					$gs_numgenes{$p[0]}++;
				}
			}
		}
		else {
			if($f3 ne '') {
				if(exists $geneset{$p[0]}) {
					push(@{$gs_genes{$p[0]}}, $p[1]);
					$gs_numgenes{$p[0]}++;
				}
			}
			else {
				push(@{$gs_genes{$p[0]}}, $p[1]);
				$gs_numgenes{$p[0]}++;
			}
		}
	}
}

print "\n";
$min_numgenes=5; $totnumsets=0; $numsetswithenoughgenes=0;
foreach(keys %gs_genes) {
	$totnumsets++;
	if($gs_numgenes{$_}>=$min_numgenes) {
		$numsetswithenoughgenes++;

		$fo=''; $fo=$_.'_'.$geneset{$_}.'.txt';
		if($fo=~/^GO:/) { $fo=~s/^GO:0*//g; }
		print "$numsetswithenoughgenes\t$fo\n";

		open HH, ">$fo";
		for($i=0;$i<=$#{$gs_genes{$_}};$i++) { print HH "${$gs_genes{$_}}[$i]\n"; }
		close HH;
	}
}

print "\nTotal no. genesets: $totnumsets\nNo. genesets with enough genes: $numsetswithenoughgenes\n\n";

close HH;

