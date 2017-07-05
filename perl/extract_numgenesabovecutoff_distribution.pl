$f1=$ARGV[0];		# Network in DAT format with the scores in descendng order
#$f2=$ARGV[1];
chomp($f3=$ARGV[1]);	# Outfile: <Score> <No.Genes> <No.Edges>

open FH,"$f1";
chomp(@f=<FH>);
close FH;
#open GH,"$f2";
#chomp(@g=<GH>);
#close GH;

open HH,">$f3";

%genes=(); $cutoff=99;
$numgenes=0; $numedges=0;
for($i=0;$i<=$#f;$i++)
{
	$_=~s///g;

	@p=();
	@p=split("\t",$f[$i]);

	if($p[2]>($cutoff/100))
	{
		$numedges++;

		unless(exists($genes{$p[0]}))
		{
			$numgenes++;
			$genes{$p[0]}++;
		}
		unless(exists($genes{$p[1]}))
		{
			$numgenes++;
			$genes{$p[1]}++;
		}
	}
	else
	{
		print HH "$cutoff\t$numgenes\t$numedges\n";
		$cutoff=$cutoff-1;
		$i--;
	}

	if(($numedges%500000) eq 0) { print "$p[0]\t$p[1]\t$p[2]\n\t$cutoff\t$numedges\t$numgenes\n"; }
}
print HH "$cutoff\t$numgenes\t$numedges\n";

close HH;
