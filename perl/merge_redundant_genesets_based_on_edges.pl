#########################################################################
# This script is designed to merge genesets that induce highly similar	#
# sets of edges in an underlying network.				#
#########################################################################

$f1=$ARGV[0];		# *.genesets
$f2=$ARGV[1];		# Edges to include; Basically a network in DAT format
$maxjac=$ARGV[2];	# Jaccard index cutoff
$mindiff=$ARGV[3];	# Set difference cutoff
chomp($f3=$ARGV[4]);	# Outfile: Redundant genesets merged

open FH,"$f1";
chomp(@gset=<FH>);
close FH;
open GH,"$f2";
chomp(@net=<GH>);
close GH;

open HH,">$f3";

#%edgescore=(); @edgearray=();
%backgroundgenelist=(); %backgroundedgelist=();
foreach (@net)
{
	@p=(); @p=split("\t",$_);

	$edge = join '__', sort($p[0], $p[1]);

	$backgroundgenelist{$p[0]}++; $backgroundgenelist{$p[1]}++;
	$backgroundedgelist{$edge}++;
}

%genegs=(); %gsnumgenes=();
foreach (@gset)
{
	@p=(); @p=split("\t",$_);

	if(exists $backgroundgenelist{$p[0]})
	{
		for($j=1;$j<=$#p;$j++)
		{
			push(@{$genegs{$p[0]}}, $p[$j]);
			$gsnumgenes{$p[$j]}++;
		}
	}
}

%gsedges=(); %gssize=(); # %gspairdesc=();
$totnumedges=0; $numedgeswithann=0; %coannedges=();
foreach $edge (keys %backgroundedgelist)
{
	$totnumedges++;
	@p=(); @p=split("__",$edge);

	if((exists $genegs{$p[0]})&&(exists $genegs{$p[1]}))
	{
		$numedgeswithann++;
		for($i=0;$i<=$#{$genegs{$p[0]}};$i++)
		{
			for($j=0;$j<=$#{$genegs{$p[1]}};$j++)
			{
				if(${$genegs{$p[0]}}[$i] eq ${$genegs{$p[1]}}[$j])
				{
					$coannedges{$edge}++;
					
					$gsid=${$genegs{$p[0]}}[$i];
					push(@{$gsedges{$gsid}}, $edge);
					$gssize{$gsid}++;
				}
			}
		}
	}
}

@gsarray = ();
@gsarray = sort {$gssize{$b} <=> $gssize{$a}} keys %gssize;

$totnumgs=0; $totnumgs = scalar(@gsarray);
$numcoannedges=0; $numcoannedges=scalar(keys %coannedges);
print "\nTot. no. genesets: $totnumgs\n\nTot. no. edges: $totnumedges\nNo. of annotated edges: $numedgeswithann\nTot. no. coannotated edges: $numcoannedges\n\n";

%redundantgs=(); %consideredgs=();
for($i=0;$i<$#gsarray;$i++)
{
	unless(exists $consideredgs{$gsarray[$i]})
	{
		$uniongstag=$gsarray[$i]; %uniongs=(); $num_red=0;
		
		%temp_edges=();
		foreach(@{$gsedges{$gsarray[$i]}}) { $temp_edges{$_}++; $uniongs{$_}++; }

		$consideredgs{$gsarray[$i]}++;

		$j=$i+1;
		while(($j <= $#gsarray)&&(($gssize{$gsarray[$i]}-$gssize{$gsarray[$j]})<$mindiff))
		{
			unless($consideredgs{$gsarray[$j]})
			{
				$common=0; $union=$gssize{$gsarray[$i]}; %temp_uniongs=%uniongs;

				foreach $edge (@{$gsedges{$gsarray[$j]}})
				{
					$temp_uniongs{$edge}++;

					if(exists $temp_edges{$edge}) { $common++; }
					else { $union++; }
				}

				$jaccard=($common/$union);
				if($jaccard >= $maxjac)
				{
					$num_red++;
					$uniongstag=$uniongstag.'__'.$gsarray[$j];
					%uniongs=%temp_uniongs;

					$consideredgs{$gsarray[$j]}++;
				}
			}
			
			$j++;
		}

		if($num_red>0)
		{
			@tags=();
			@tags=split("__",$uniongstag);

			foreach $gs (@tags)
			{
				delete $gsedges{$gs};
			}
			
			$num_edges_in_union=0;
			foreach $edge (keys %uniongs)
			{
				$num_edges_in_union++;
				push(@{$gsedges{$uniongstag}}, $edge);
			}
		}
	}
}

%remainingedges=(); %remaininggs=();
foreach $gs (keys %gsedges)
{
	if($gs=~/__/)
	{
		@p=();
		@p=split("__",$gs);

		foreach $edge (@{$gsedges{$gs}})
		{
			#print HH "$edge\t$p[0]\n";
			$remainingedges{$edge}++;
			$remaininggs{$p[0]}++;
		}
	}
	else
	{
		foreach $edge (@{$gsedges{$gs}})
		{
			#print HH "$edge\t$gs\n";
			$remainingedges{$edge}++;
			$remaininggs{$gs}++;
		}
	}
}

$nummergedgs=0; $nummergedgs = scalar(keys %gsedges);
$numremainingedges=0; $numremainingedges = scalar(keys %remainingedges);
print "No. merged genesets: $nummergedgs\nNo. remaining edges: $numremainingedges\n\n";

foreach(keys %remaininggs) { print HH "$_\n"; }
close HH;

