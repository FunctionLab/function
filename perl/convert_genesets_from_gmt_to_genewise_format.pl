$f1=$ARGV[0];		# GeneSet in GMT format (typically from MSigDB); <GS> <gene1> <gene2> ...
chomp($f3=$ARGV[1]);	# Output: GeneSet in gene-wise format: <Gene> <GS1> <GS2> ...

open FH,"$f1";
chomp(@f=<FH>);
close FH;

open HH,">$f3";

%genegs=();
foreach(@f)
{
	$_=~s///g;

	@p=();
	@p=split("\t",$_);

	if(($#p eq 1)&&($p[1]=~/\,/))
	{
		@q=();
		@q=split(",",$p[1]);

		for($j=0;$j<=$#q;$j++)
		{
			push(@{$genegs{$q[$j]}}, $p[0]);
		}
	}
	else
	{
		for($j=2;$j<=$#p;$j++)
		{
			push(@{$genegs{$p[$j]}}, $p[0]);
		}
	}
}

foreach $gene (sort {$a cmp $b} keys %genegs)
{
	print HH "$gene";

	for($j=0;$j<=$#{$genegs{$gene}};$j++)
	{
		print HH "\t${$genegs{$gene}}[$j]";
	}

	print HH "\n";
}

close HH;
