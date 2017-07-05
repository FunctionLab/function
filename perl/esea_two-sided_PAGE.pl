#########################################################################################
# This script calculates functional enrichment on a edge list using Parametric Geneset	#
# Enrichment Analysis (PAGE) described in Kim & Volsky (2005) BMC Bioinformatics.	#									#
# 											#
# Also refer to http://expressome.kobic.re.kr/GAzer/supplement.jsp for pointers on the	#
# performance of this approach compared to others.					#
# 											#
# The input is a list or a matrix of all the assayed edges (first column) and their	#
# attributes or	scores in the rest of the columns. This data matrix should contain	#
# weighted-correlation-changes, signed negative log p-values or z-scores of all the	#
# assayed edges, the column containing which can be specified. Since these values can	#
# be ranked to produce a two-tailed distribution, enrichment will also be a two-sided	#
# test to check if a given edgeset is has a significantly high or low average score.	#
# 											#
# If the matrix has multiple columns of test statistics from assays/experiments, one	#
# can also provide more than one column for analysis, which are	then considered		#
# sequentially. If all the columns have to be considered, then provide no input. If	#
# some columns are to be chosen, then the input can be something like "3_5:8_11_13",	#
# which means columns 3, 5, 6, 7, 8, 11 & 13.						#
# 											#
# The output is a matrix containing the edgesets along the rows and specified columns	#
# with each cell (i,j) containing the enrichment score (-1)log10(Q-value) of edgeset i	#
# in the ranking provided by column j. Q-values come from corrections for multiple	#
# hypothesis testing using Benjamini-Hochberg method.					#
#											#
# Possible Extenstions:									#
# ====================									#
# i) Empirical distributions are build based on a permutation test (k=10,000) for	#
# edgesets of size less than 10.							#
# ii) P-values computed here are also stored for conversion to q-values using methods	#
# other than BH (implemented here) like the Storey, Tishirani (2003) PNAS method in R.	#
#########################################################################################

use Time::SoFar qw( runtime runinterval figuretimes );
use List::Util 'shuffle';
use Statistics::Distributions;

# Reading in options and all files
##################################

$elapsed=0; $elapsed = runtime();
print "\n$elapsed: Reading in options and files ...";

$f1=$ARGV[0];		# Input edges; Contains header; Can be multi-column matrix;
$coltag=$ARGV[1];	# Columns containing test statistic; Index begins at 0; if 'all', then all columns are considered
$anal=$ARGV[2];		# Analysis: 'intra' or 'inter'
$tail=$ARGV[3];		# Tail: 'twotail', 'utail' or 'ltail'
$misse=$ARGV[4];	# Treatment of missing edges: 'i' (ignore), or '0' (consider, and assign a weight of zero)
$f2=$ARGV[5];		# Functional annotations = GeneSets; *.genesets
#$f3=$ARGV[6];		# GeneSet relationships: DAG-relationships for GO: <GOID> <GOID> <Distance_in_DAG>; Could be 'NA'
chomp($f4=$ARGV[6]);	# List of genesets to include in the analysis

open FH,"$f1";
chomp(@f = grep {!/^#/} <FH>);
close FH;
open GH,"$f2";
chomp(@gset = grep {!/^#/} <GH>);
close GH;
open IH,"$f4";
chomp(@gslist = grep {!/^#/} <IH>);
close IH;

$gsetdescfile=''; $gsetdescfile=$f2.'.desc';

open DH2,"$gsetdescfile";
chomp(@gsdesc=<DH2>);
close DH2;

$fout1 = $f2; $fout1 =~ s/\.\.\///g; $fout1 =~ s/^.*\///g; $fout1 =~s/\.genesets//g;
$fout1 = $f1.'-'.$fout1.'.'.$anal.'-esea.zscore.mat';
$fout2 = $fout1; $fout2 =~ s/\.zscore\.mat$/\.esfdr\.mat/g;
open JH,">$fout1";
open HH,">$fout2";


@colarray=();
if($coltag eq 'all')
{
	@p=();
	@p=split("\t",$f[0]);

	for($i=2;$i<=$#p;$i++)
	{
		push(@colarray, $i);
	}
}
else
{
	@p=(); @p=split("_",$coltag); $k=0;

	foreach(@p)
	{
		if($_=~/:/)
		{ 
			@q=();
			@q=split(":",$_);

			for($i=$q[0];$i<=$q[1];$i++)
			{
				splice(@colarray, $k, 0, $i);
				$k++;
			} 
		} 
		else
		{ 
			push(@colarray, $_);
			$k++;
		} 
	}
}

print "\nInput file: $f1\nOutput file: $f3\nColumns: @colarray\n";


# Indexing input edges based on background lists and obtaining their conditional scores
##########################################################################################

$elapsed=0; $elapsed = runtime();
print "\n$elapsed: Indexing input edges, genesets desc and background lists ...\n";

%edgescore=(); @edgearray=(); %backgroundgenelist=(); %backgroundedgelist=();

@sum_edgescore=(); # @sqsum_edgescore=();
foreach(@colarray) { push(@sum_edgescore, 0); } # push(@sqsum_edgescore, 0); }

for($i=1;$i<=$#f;$i++)
{
	$f[$i]=~s///g;

	@p=();
	@p=split("\t",$f[$i]);

	$edge = join '__', sort($p[0], $p[1]);

	push(@edgearray, $edge);
	$backgroundgenelist{$p[0]}++; $backgroundgenelist{$p[1]}++;
	$backgroundedgelist{$edge}++;

	for($j=0;$j<=$#colarray;$j++)
	{
		push(@{$edgescore{$edge}}, $p[$colarray[$j]]);
		$sum_edgescore[$j]+=$p[$colarray[$j]];
	}
}

$univ=0; $univ=scalar(@edgearray);

if($misse eq '0')
{
	$num_missingedges=0;
	@all_genes = keys %backgroundgenelist;

	for($i=0;$i<$#all_genes;$i++)
	{
		for($j=($i+1);$j<=$#all_genes;$j++)
		{
			$edge = join '__', sort($all_genes[$i], $all_genes[$j]);

			unless(exists $backgroundedgelist{$edge})
			{
				$num_missingedges++;
				$univ++;
				$backgroundedgelist{$edge}++;
			}
		}
	}
}

@mean_edgescore=();
for($j=0;$j<=$#colarray;$j++) { push(@mean_edgescore, ($sum_edgescore[$j]/$univ)); }

print "\n\tUniv: $univ";
print "\n\tSum of edgescores: @sum_edgescore"; #\n\tUniv: $univ";

@sd_edgescore=();
for($j=0;$j<=$#colarray;$j++)
{
	$sqsum=0;
	foreach(@edgearray)
	{
		$sqsum+=(${$edgescore{$_}}[$j]-$mean_edgescore[$j])**2;
	}
	if($misse eq '0') { $sqsum+=($num_missingedges*($mean_edgescore[$j]**2)); }

	push(@sd_edgescore, sqrt($sqsum/$univ));
}

print "\n\tMean edgescore: @mean_edgescore\n\tSD edgescores: @sd_edgescore\n";


# Indexing gs descriptions and member edges
################################################

$elapsed=0; $elapsed = runtime();
print "\n$elapsed: Indexing gs descriptions and edges ...";

%gsdesc=();
foreach(@gsdesc)
{
	$_=~s///g;

	@p=();
	@p=split("\t",$_);

	$gsdesc{$p[0]}=$p[1];
}

%gs_toincl=();
foreach(@gslist) { $gs_toincl{$_}++; }

%gsgenes=(); %genegs=(); %gsnumgenes=();
foreach $gene (@gset)
{
	$gene=~s///g;

	@p=();
	@p=split("\t",$gene);

	if(exists $backgroundgenelist{$p[0]})
	{
		for($j=1;$j<=$#p;$j++)
		{
			if(exists $gs_toincl{$p[$j]})
			{
				push(@{$genegs{$p[0]}}, $p[$j]);
				push(@{$gsgenes{$p[$j]}}, $p[0]);
				$gsnumgenes{$p[$j]}++;
			}
		}
	}
}

%gsedges=(); %gssize=(); %gspairdesc=(); %coannotated_edges=();
foreach $gs (keys %gsgenes)
{
	$gspairdesc{$gsid}=$gsdesc{$gsid};
	
	for($i=0;$i<$#{$gsgenes{$gs}};$i++)
	{
		for($j=($i+1);$j<=$#{$gsgenes{$gs}};$j++)
		{
			$edge = join '__', sort(${$gsgenes{$gs}}[$i], ${$gsgenes{$gs}}[$j]);
			
			if(exists $backgroundedgelist{$edge})
			{
				push(@{$coannotated_edges{$edge}}, $gs);

				if($anal eq 'intra')
				{
					push(@{$gsedges{$gsid}}, $edge);
					$gssize{$gsid}++;
				}
			}
		}
	}
}

if($anal eq 'inter')
{
	foreach $edge (keys %backgroundedgelist)
	{
		@p=();
		@p=split("__",$edge);

		if((not exists $coannotated_edges{$edge})&&(exists $genegs{$p[0]})&&(exists $genegs{$p[1]}))
		{
			for($k=0;$k<=$#{$genegs{$p[0]}};$k++)
			{
				for($l=0;$l<=$#{$genegs{$p[1]}};$l++)
				{
					@gsorder = sort(${$genegs{$p[0]}}[$k], ${$genegs{$p[1]}}[$l]);
					$gsid = join '__', @gsorder;
					push(@{$gsedges{$gsid}}, $edge);

					$gssize{$gsid}++;
					unless(exists $gspairdesc{$gsid})
					{
						$gspairdesc{$gsid}=$gsdesc{$gsorder[0]}.'|'.$gsdesc{$gsorder[1]};
					}
				}
			}
		}
	}
}

@gsarray=(); $totnumgs=0; $numgswithenoughedges=0;
$minnumedges=10;

$totnumgs=scalar(keys %gssize);

foreach $gs (sort {$gssize{$b}<=>$gssize{$a}} keys %gssize)
{
	if($gssize{$gs}>=$minnumedges)
	{
		$numgswithenoughedges++;
		push(@gsarray, $gs);

		if($anal eq 'inter')
		{
			%gsgenes=();
			foreach $e (@{$gsedges{$gs}})
			{
				@p=(); @p=split("__",$e);
				$gsgenes{$p[0]}++; $gsgenes{$p[1]}++;
			}
			$gsnumgenes{$gs}=scalar(keys %gsgenes);
		}
	}
	else
	{
		last;
	}
}

print "\n\tTotal no. of gss: $totnumgs\n\tNo. of gss with $minnumedges or more edges: $numgswithenoughedges\n\n";


# Calculating edge set scores, enrichment z-scores and their pvalues
####################################################################

print JH "GS.ID\tDesc\tNo.Genes\tNo.Edges\tSum.Score\tMean.Score\tZ.Score\n";

#@p=();
#@p=split("\t",$f[0]);

#for($j=0;$j<=$#colarray;$j++)
#{
#	print JH "\t$p[$colarray[$j]]";
#}

#print JH "\n";

$elapsed=0; $elapsed = runtime();
print "$elapsed: Calculating gs scores ...", "\n";
# print " ...\n";

%gsscore=();
# %rand_gsscores=(); %emp_mean=(); %emp_sd=(); $n=10000;
# %gsparzscore=(); %par_sd=();
%gsparpval=(); $numtests=0;

for($i=0;$i<=$#gsarray;$i++)
{
	$numtests++;
	$gssize = $gssize{$gsarray[$i]};

	print JH "$gsarray[$i]\t$gspairdesc{$gsarray[$i]}\t$gsnumgenes{$gsarray[$i]}\t$gssize";

	# Mean edgeset score
	@{$gsscore{$gsarray[$i]}} = @{meangsscore(\@{$gsedges{$gsarray[$i]}}, $gssize)};
	
	@{$gsparpval{$gsarray[$i]}}=();		# P-value

	for($j=0;$j<=$#colarray;$j++)
	{
		$meangsscore = ${$gsscore{$gsarray[$i]}}[$j];
		$sumgsscore = sprintf("%.3g", ($meangsscore*$gssize));
		
		$parsd = ($sd_edgescore[$j]/sqrt($gssize));

		$zscore = sprintf("%.3f", ($meangsscore - $mean_edgescore[$j])/$parsd);
		# push(@{$gsparzscore{$gsarray[$i]}}, $zscore);

		if($tail eq 'twotail') { $pvalue = sprintf("%.5g", (Statistics::Distributions::uprob (abs($zscore)))*2); }
		elsif($tail eq 'utail') { $pvalue = sprintf("%.5g", (Statistics::Distributions::uprob ($zscore))); }
		elsif($tail eq 'ltail') { $pvalue = sprintf("%.5g", (1 - Statistics::Distributions::uprob ($zscore))); }
		push(@{$gsparpval{$gsarray[$i]}}, $pvalue);
	
		$meantoprint = sprintf("%.3f", $meangsscore);
		# print JH "\t$meantoprint\t$mean_edgescore[$j]\t$parsd\t$zscore\t$pvalue";
		# print JH "\t$sumgsscore\t$meangsscore\t$zscore";
		print JH "\t$sumgsscore\t$meantoprint\t$zscore";
	}

	print JH "\n";
}


# Multiple Hypotheses Testing correction
############################################################

print HH "GS.ID\tDesc\tNo.Genes\tNo.Edges\tSum.Score\tMean.Score\tES.FDR\n";

@p=();
@p=split("\t",$f[0]);

$elapsed=0; $elapsed = runtime();
print "\n$elapsed: Calculating corrected qvalues ...\n";

%gsesparfdr=();
for($j=0;$j<=$#colarray;$j++)
{
	#print HH "\t$p[$colarray[$j]]";
	
	$minparpval=1;
	foreach $gs (keys %gsparpval)
	{
		if((${$gsparpval{$gs}}[$j] ne 0)&&(${$gsparpval{$gs}}[$j] < $minparpval))
		{
			$minparpval=${$gsparpval{$gs}}[$j];
		}
	}

	$rank=$numtests;
	foreach (sort {${$gsparpval{$b}}[$j] <=> ${$gsparpval{$a}}[$j]} keys %gsparpval)
	{
		if(${$gsparpval{$_}}[$j] eq 0) { ${$gsparpval{$_}}[$j] = 0.1*sprintf("%.6g", $minparpval); }

		$fdr = ${$gsparpval{$_}}[$j]*($numtests/$rank); if($fdr>1) { $fdr = 1; }
		
		$sign=1; if(${$gsscore{$_}}[$j]<0) { $sign=-1; }
		$esparfdr = $sign*(-1)*log($fdr)/log(10);
		
		push(@{$gsesparfdr{$_}}, $esparfdr);

		$rank--;
	}
}

#print HH "\n";

# Printing results
##################

$elapsed=0; $elapsed = runtime();
print "$elapsed: Printing results ...\n";

foreach $gs (keys %gsesparfdr)
{
	print HH "$gs\t$gspairdesc{$gs}\t$gsnumgenes{$gs}\t$gssize{$gs}\t\t";

	for($j=0;$j<=$#colarray;$j++)
	{
		$esfdrtoprint = sprintf("%.3f", ${$gsesparfdr{$gs}}[$j]);
		print HH "\t$esfdrtoprint";
	}

	print HH "\n";
}

$elapsed=0; $elapsed = runtime();
print "$elapsed: DONE ...\n\n";

close HH;


# Subroutine to calculate the mean score of a gs
#####################################################

sub meangsscore
{
	my ($tempga, $size) = @_;
	
	my @totscore=(); foreach(@colarray) { push(@totscore, 0); }
	
	foreach (@$tempga)
	{
		for(my $k=0;$k<=$#colarray;$k++)
		{
			my $score=0;
			if(exists $edgescore{$_}) { $score=${$edgescore{$_}}[$k]; }
			$totscore[$k]+=$score;
		}
	}

	my @meanscore=();
	for($k=0;$k<=$#colarray;$k++) { push(@meanscore, ($totscore[$k]/$size)); }

	return \@meanscore;
}
