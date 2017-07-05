#########################################################################
# This script takes z-scores reported by ESEA, converting them to       #
# p-values, which will be used as weights to obtain a set-cover	of all	#
# the annotated edges in a network.                                     #
# 									                                    #
# If provided with only intra-set z-scores, the set-cover will be       #
# calculated over the set of coannotated edges. However, if inter-set	#
# z-scores are also provided, all annotated edges will be taken into	#
# account.								                                #
#########################################################################

use Time::SoFar qw( runtime runinterval figuretimes );
use Statistics::Distributions;

$f1=$ARGV[0];		# ESEA ouput: <geneset> ... <z-score>
$f2=$ARGV[1];		# Geneset annotations: .genesets file
$anal=$ARGV[2];		# 'intra' or 'both'
$tail=$ARGV[3];		# tail: 'twotail', 'utail', or 'ltail'
$f3=$ARGV[4];		# Network in DAT format: <nodeA> <nodeB> <weight>
chomp($f4=$ARGV[5]);	# Output: Collection of sets and their scores

$elapsed = runtime();
print "\n$elapsed: Reading in files ... ";

open GS,"$f1";
chomp(@gs=<GS>);
close GS;
open GG,"$f2";
chomp(@gg=<GG>);
close GG;
open NH,"$f3";
chomp(@net=<NH>);
close NH;

open HH,">$f4";

$elapsed = runtime();
print "Done\n$elapsed: Recording background, and geneset scores, genes ... ";

%backgroundgenes=(); %backgroundedges=();
foreach(@net)
{
	@p=(); @p=split("\t",$_);

	$backgroundgenes{$p[0]}++; $backgroundgenes{$p[1]}++;
	$edge = join '__', sort($p[0], $p[1]);
	$backgroundedges{$edge}++;
}

print HH "$gs[0]\tNo.New.Edges\tFrac.Edges.Covered\tCostBenf.Ratio\n";
%gs_score=(); %gs_info=();
for($i=1;$i<=$#gs;$i++)
{
	@p=();
	@p=split("\t",$gs[$i]);

	if($tail eq 'twotail') { $pval = sprintf("%.6g", (Statistics::Distributions::uprob (abs($p[$#p])))*2); }
	elsif($tail eq 'utail') { $pval = sprintf("%.6g", Statistics::Distributions::uprob ($p[$#p])); }
	elsif($tail eq 'ltail') { $pval = sprintf("%.6g", (1 - Statistics::Distributions::uprob ($p[$#p]))); }

	$gs_score{$p[0]}=$pval;
	$gs_info{$p[0]}=$gs[$i];
}

%genegs=(); %gsgenes=();
foreach(@gg)
{
	@p=();
	@p=split("\t",$_);

	if(exists $backgroundgenes{$p[0]})
	{
		for($j=1;$j<=$#p;$j++)
		{
			push(@{$genegs{$p[0]}}, $p[$j]);
			push(@{$gsgenes{$p[$j]}}, $p[0]);
		}
	}
}

$elapsed = runtime();
print "Done\n$elapsed: Populating edge sets ... ";

foreach $edge (keys %backgroundedges)
{
	@p=();
	@p=split("__",$edge);

	if((exists $genegs{$p[0]})&&(exists $genegs{$p[1]}))
	{
		for($i=0;$i<=$#{$genegs{$p[0]}};$i++)
		{
			for($j=0;$j<=$#{$genegs{$p[1]}};$j++)
			{
				if(${$genegs{$p[0]}}[$i] eq ${$genegs{$p[1]}}[$j])
				{
					$gsid=${$genegs{$p[0]}}[$i];
					push(@{$coannotated_edges{$edge}}, $gsid);

					push(@{$gsedges{$gsid}}, $edge);
					$gssize{$gsid}++;
				}
			}
		}

		if(($anal eq 'both')&&(not exists $coannotated_edges{$edge}))
		{
			for($i=0;$i<=$#{$genegs{$p[0]}};$i++)
			{
				for($j=0;$j<=$#{$genegs{$p[1]}};$j++)
				{
					$gsid = join '|', sort(${$genegs{$p[0]}}[$i], ${$genegs{$p[1]}}[$j]);

					push(@{$gsedges{$gsid}}, $edge);
					$gssize{$gsid}++;
				}
			}
		}
	}
}

%ann_edges=();
foreach $gs (keys %gsedges)
{
	if($gssize{$gs}<10) { delete $gsedges{$gs}; }
	else { foreach $edge (@{$gsedges{$gs}}) { $ann_edges{$edge}++; }  }
}

$elapsed = runtime();
print "Done\n$elapsed: Performing set-cover ...\n";

$tot_num_edges=scalar(keys %ann_edges);
%uncovered_edges=%ann_edges; $num_uncovered_edges=scalar(keys %ann_edges);
%unselected_sets=%gsedges; %selected_sets=(); $num_covered_edges=0;

print "Tot. no. uncovered edges: $num_uncovered_edges\n"; $count=0;

while(($num_uncovered_edges > 0)&&($count<3))
{
	$count++;
	$mincostbenf=1; $num_new_edges=0; $innercount=0;
	foreach $gs (keys %unselected_sets)
	{
		$innercount++;
		$tempnum_new_edges=0;
		foreach $edge (@{$gsedges{$gs}})
		{
			if(exists $uncovered_edges{$edge})
			{
				$tempnum_new_edges++;
			}
		}

		if($tempnum_new_edges>0)
		{
			$ratio=($gs_score{$gs}/$tempnum_new_edges);
			if($innercount eq 1) { if(exists $gs_score{$gs}) { print "\t$gs\t$gs_score{$gs}\t$tempnum_new_edges\t$ratio\n"; } else { print "\t$gs\tNo Score !!\n"; } }

			if($ratio < $mincostbenf) { $mincostbenf=$ratio; $mings=$gs; $num_new_edges=$tempnum_new_edges; }
		}
	}

	$num_covered_edges+=$num_new_edges;
	$mincostbenf_toprint=sprintf("%.5g", $mincostbenf);
	$frac_edges_covered=sprintf("%.5g", ($num_covered_edges/$tot_num_edges));

	print HH "$gs_info{$mings}\t$num_new_edges\t$frac_edges_covered\t$mincostbenf_toprint\n";

	$selected_sets{$mings}=$mincostbenf;
	delete $unselected_sets{$mings};
	foreach $edge (@{$gsedges{$mings}}) { if(exists $uncovered_edges{$edge}) { delete $uncovered_edges{$edge}; } }
	$num_uncovered_edges=scalar(keys %uncovered_edges);

	print "\t$gs_info{$mings}\n";
	print "\t\tNo. new edges added: $num_new_edges\tCost-Benefit ratio: $mincostbenf_toprint\n";
	print "\t\tFrac. edges covered: $frac_edges_covered\tNo. uncovered edges: $num_uncovered_edges\n";
}

$elapsed = runtime();
print "$elapsed: DONE\n\n";

close HH;
