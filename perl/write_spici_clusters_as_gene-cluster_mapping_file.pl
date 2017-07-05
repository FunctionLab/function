$f1=$ARGV[0];		# .clusters file; spici output
$f2=$ARGV[1];		# Total list of genes
$par=$ARGV[2];		# Processing parameter: 'g' for cluster-genesets format OR 'n' for cluster-based network format OR 'c' for cluster-based context files
chomp($f3=$ARGV[3]);	# Oufile: Gene-ClusterID OR network in DAT format with edge weights equal to 1 OR directory to output contexts

open FH,"$f1";
chomp(@f=<FH>);
close FH;
open GH,"$f2";
chomp(@g=<GH>);
close GH;

if(($par eq 'g')||($par eq 'n')) { open HH,">$f3"; }

%clusteredgenes=(); $numedges=0;
for($i=0;$i<=$#f;$i++)
{
	$f[$i]=~s///g;

	@p=();
	@p=split("\t",$f[$i]);

	if($par eq 'g')
	{
		$clusternum=sprintf("%04d", ($i+1)); $clusterid='C'.$clusternum;

		for($j=0;$j<=$#p;$j++)
		{
			$clusteredgenes{$p[$j]}++;
			print HH "$clusterid\t$p[$j]\n";
		}
	}
	elsif($par eq 'n')
	{
		for($j=0;$j<$#p;$j++)
		{
			for($k=$j+1;$k<=$#p;$k++)
			{
				$numedges++;
				print HH "$p[$j]\t$p[$k]\t1\n";
			}
		}
	}
	elsif($par eq 'c')
	{
		$clusternum=sprintf("%04d", ($i+1)); $clusterfn=$f3.'C'.$clusternum;
		print "$clusterfn ...\n";
		open HH, ">$clusterfn";

		for($j=0;$j<=$#p;$j++)
		{
			#$clusteredgenes{$p[$j]}++;
			print HH "$p[$j]\n";
		}

		close HH;
	}
}

#if($par eq 'g')
#{
#	foreach(@g)
#	{
#		$_=~s///g;
#
#		@p=();
#		@p=split("\t",$_);
#
#		unless(exists $clusteredgenes{$p[0]})
#		{
#			print HH "$p[0]\tCluster0000\n";
#		}
#	}
#}

if($par eq 'n') { print "\nNo. of edges: $numedges\n\n"; }
elsif(($par eq 'g')||($par eq 'n')) { close HH; }
