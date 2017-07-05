$f1=$ARGV[0];		# File containing list of DChecker outfiles
#$f2=$ARGV[1];
chomp($f3=$ARGV[1]);	# Outfile containing a single DChecker-like output

open FH,"$f1";
chomp(@f=<FH>);
close FH;
#open GH,"$f2";
#chomp(@g=<GH>);
#close GH;

open HH,">$f3";

$numpos=0; $numneg=0;
%cut_info=();
foreach(@f)
{
	@d=();
	open DH, "$_";
	chomp(@d=<DH>);
	close DH;

	print "$_ $#d ...\n";

	@p=(); @p=split("\t",$d[0]); $numpos+=$p[2];
	@p=(); @p=split("\t",$d[1]); $numneg+=$p[2];
	$header=$d[2];

	for($i=3;$i<$#d;$i++)
	{
		@p=();
		@p=split("\t",$d[$i]);
		print "$i\t@p\n";
		
		$par=0;
		if(exists $cut_info{$p[0]})
		{
			for($j=1;$j<=$#p;$j++) { ${$cut_info{$p[0]}}[$j-1]+=$p[$j]; }
		}
		else
		{
			print "\t$#p $p[0]"; $par=1;
			for($j=1;$j<=$#p;$j++) { push(@{$cut_info{$p[0]}}, $p[$j]); print " $j:$p[$j]"; }
		}
		if($par=1) { print "\n"; }
	}
}

print HH "#\tP\t$numpos\n#\tN\t$numneg\n$header\n";
foreach $cut (sort {$a <=> $b} keys %cut_info)
{
	print HH "$cut\t$cut_info{$cut}\n";
}

close HH;
