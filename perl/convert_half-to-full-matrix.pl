$f1=$ARGV[0];		# Half-matrix
chomp($f3=$ARGV[1]);	# Outfile: Full-matrix

open FH,"$f1";
chomp(@f=<FH>);
close FH;

open HH,">$f3";

$f[0]=~s///g;
@p=();
@p=split("\t",$f[0]);

%colidx=(); $numcol=0; @colarray=();
for($j=1;$j<=$#p;$j++)
{
	$colidx{$p[$j]}=$j-1;
	$numcol++;

	push(@colarray, $p[$j]);
	print HH "\t$p[$j]";
}
print HH "\n";

print "\nNo. of columns: $numcol\n\n";

for($i=1;$i<=$#f;$i++)
{
	$f[$i]=~s///g;

	@p=();
	@p=split("\t",$f[$i]);

	for($j=1;$j<=$#p;$j++)
	{
		$mat[$colidx{$p[0]}][$j-1]=$p[$j];
	}
}

for($i=0;$i<=$#colarray;$i++)
{
	for($j=0;$j<=$#colarray;$j++)
	{
		if($mat[$i][$j] eq '')
		{
			$mat[$i][$j]=$mat[$j][$i];
		}
	}
}

for($i=0;$i<=$#colarray;$i++)
{
	print HH "$colarray[$i]";
	for($j=0;$j<=$#colarray;$j++)
	{
		print HH "\t$mat[$i][$j]";
	}
	print HH "\n";
}

close HH;
