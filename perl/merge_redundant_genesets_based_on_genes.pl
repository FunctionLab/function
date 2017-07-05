$f1=$ARGV[0];		# *.genesets
$f2=$ARGV[1];		# Genes to include
$maxjac=$ARGV[2];	# Jaccard index cutoff
$mindiff=$ARGV[3];	# Set difference cutoff
chomp($f3=$ARGV[4]);	# Outfile: Redundant genesets merged

open FH,"$f1";
chomp(@f=<FH>);
close FH;
open GH,"$f2";
chomp(@g=<GH>);
close GH;

open HH,">$f3";

%genes_to_incl=();
foreach(@g) { $genes_to_incl{$_}++; }

%gs_genes=(); %gs_size=(); %all_genes=();
foreach(@f)
{
	@p=();
	@p=split("\t",$_);

	if(exists $genes_to_incl{$p[0]})
	{
		for($j=1;$j<=$#p;$j++)
		{
			$gs_size{$p[$j]}++;
			push(@{$gs_genes{$p[$j]}}, $p[0]);

			$all_genes{$p[0]}++;
		}
	}
}

@gs_array = ();
@gs_array = sort {$gs_size{$b} <=> $gs_size{$a}} keys %gs_size;

$totnumgs=0; $totnumgs = scalar(@gs_array);
$totnumgenes=0; $totnumgenes = scalar(keys %all_genes);
print "\nTot. no. genesets: $totnumgs\nTot. no. genes: $totnumgenes\n\n";

%redundant_gs=(); %considered_gs=(); # $mindiff=3;
for($i=0;$i<$#gs_array;$i++)
{
	unless(exists $considered_gs{$gs_array[$i]})
	{
		$union_gstag=$gs_array[$i]; %union_gs=(); $num_red=0;
		
		%temp_genes=();
		foreach(@{$gs_genes{$gs_array[$i]}}) { $temp_genes{$_}++; $union_gs{$_}++; }

		$considered_gs{$gs_array[$i]}++;

		$j=$i+1;
		while(($j <= $#gs_array)&&(($gs_size{$gs_array[$i]}-$gs_size{$gs_array[$j]})<$mindiff))
		{
			unless($considered_gs{$gs_array[$j]})
			{
				$common=0; $union=$gs_size{$gs_array[$i]}; %temp_union_gs=%union_gs;

				foreach $gene (@{$gs_genes{$gs_array[$j]}})
				{
					$temp_union_gs{$gene}++;

					if(exists $temp_genes{$gene}) { $common++; }
					else { $union++; }
				}

				$jaccard=($common/$union);
				if($jaccard >= $maxjac)
				{
					$num_red++;
					$union_gstag=$union_gstag.'__'.$gs_array[$j];
					%union_gs=%temp_union_gs;

					$considered_gs{$gs_array[$j]}++;
				}
			}
			
			$j++;
		}

		if($num_red>0)
		{
			@tags=();
			@tags=split("__",$union_gstag);

			foreach $gs (@tags)
			{
				delete $gs_genes{$gs};
			}
			
			$num_genes_in_union=0;
			foreach $gene (keys %union_gs)
			{
				$num_genes_in_union++;
				push(@{$gs_genes{$union_gstag}}, $gene);
			}
		}
	}
}

%remaining_genes=();
foreach $gs (keys %gs_genes)
{
	if($gs=~/__/)
	{
		@p=();
		@p=split("__",$gs);

		foreach $gene (@{$gs_genes{$gs}})
		{
			print HH "$gene\t$p[0]\n";
			$remaining_genes{$gene}++;
		}
	}
	else
	{
		foreach $gene (@{$gs_genes{$gs}})
		{
			print HH "$gene\t$gs\n";
			$remaining_genes{$gene}++;
		}
	}
}

$nummergedgs=0; $nummergedgs = scalar(keys %gs_genes);
$numremaininggenes=0; $numremaininggenes = scalar(keys %remaining_genes);
print "No. merged genesets: $nummergedgs\nNo. remaining genes: $numremaininggenes\n\n";

close HH;

