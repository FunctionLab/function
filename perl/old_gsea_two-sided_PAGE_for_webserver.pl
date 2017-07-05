#########################################################################################
# This script calculates functional enrichment on a gene list using Parametric Geneset  #
# Enrichment Analysis (PAGE) described in Kim & Volsky (2005) BMC Bioinformatics.       #
#                                                                                       #
# Also refer to http://expressome.kobic.re.kr/GAzer/supplement.jsp for pointers on the  #
# performance of this approach compared to others.                                      #
#                                                                                       #
# The input is a list or a matrix of all the assayed genes (first column) and their     #
# attributes or	scores in the rest of the columns. This data matrix should contain      #
# t-statistics, fold-changes, signed negative log p-values or z-scores of all the       #
# assayed genes, the column containing which can be specified. Since these values can   #
# be ranked to produce a two-tailed distribution, enrichment will also be a two-sided   #
# test to check if a given geneset is has a significantly high or low average score.    #
#                                                                                       #
# If the matrix has multiple columns of test statistics from assays/experiments, one    #
# can also provide more than one column for analysis, which are	then considered         #
# sequentially. If all the columns have to be considered, then provide no input. If     #
# some columns are to be chosen, then the input can be something like "3_5:8_11_13",    #
# which means columns 3, 5, 6, 7, 8, 11 & 13.                                           #
#                                                                                       #
# The output is a matrix containing the genesets along the rows and specified columns   #
# with each cell (i,j) containing the enrichment score (-1)log10(Q-value) of geneset i  #
# in the ranking provided by column j. Q-values come from corrections for multiple      #
# hypothesis testing using Benjamini-Hochberg method.                                   #
#                                                                                       #
# Possible Extenstions:                                                                 #
# ====================                                                                  #
# i) Empirical distributions are build based on a permutation test (k=10,000) for       #
# genesets of size less than 10.                                                        #
# ii) P-values computed here are also stored for conversion to q-values using methods   #
# other than BH (implemented here) like the Storey, Tishirani (2003) PNAS method in R.  #
#########################################################################################

use Time::SoFar qw( runtime runinterval figuretimes );
use List::Util 'shuffle';
use Statistics::Distributions;
use strict;

# Reading in options and all files
##################################

my $elapsed=0; $elapsed = runtime();
print "\n$elapsed: Reading in options and files ...";

my($f1, $coltag, $f2, $qvalcutoff)=@ARGV;
#$f1:		Input genes; Contains header; Can be multi-column matrix;
#$coltag:	Columns containing test statistic; Index begins at 0; if 'all', then all columns are considered
#$f2:		Functional annotations = GeneSets; *.genesets
#$qvalcutoff:	Qvalue cutoff used for deciding significant genesets

my(@f, @gset)
open FH,"$f1"; chomp(@f=<FH>); close FH;
open GH,"$f2"; chomp(@gset=<GH>); close GH;

my $gsetdescfile=''; $gsetdescfile=$f2.'.desc';
open DH2,"$gsetdescfile"; chomp(my @gsdesc=<DH2>); close DH2;

my @colarray=(); my (@p, @q);
if($coltag eq 'all')
{
	@p=(); @p=split("\t",$f[0]);
	for(my $i=1;$i<=$#p;$i++) { push(@colarray, $i); }
}
else
{
	@p=(); @p=split("_",$coltag); my $k=0;
	foreach(@p)
	{
		if($_=~/:/)
		{ 
			@q=(); @q=split(":",$_);
			for(my $i=$q[0];$i<=$q[1];$i++)
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

print "\nInput file: $f1\nOutput file: $f3\n";


# Indexing input genes based on background lists and obtaining their conditional scores
##########################################################################################

$elapsed=0; $elapsed = runtime();
print "\n$elapsed: Indexing input genes, genesets desc and background lists ...";

my %genescore=(); my @genearray=(); my %backgroundlist=();

my @sum_genescore=(); # @sqsum_genescore=();
foreach(@colarray) { push(@sum_genescore, 0); } # push(@sqsum_genescore, 0); }

#print "\n\tSum of genescore: @sum_genescore";
for(my $i=1;$i<=$#f;$i++)
{
	$f[$i]=~s///g;
	@p=(); @p=split("\t",$f[$i]);

	push(@genearray, $p[0]);
	$backgroundlist{$p[0]}++;

	for(my $j=0;$j<=$#colarray;$j++)
	{
		push(@{$genescore{$p[0]}}, $p[$colarray[$j]]);

		$sum_genescore[$j]+=$p[$colarray[$j]];
	}
}

my $univ=0; $univ=scalar(@genearray);

my @mean_genescore=();
for(my $j=0;$j<=$#colarray;$j++) { push(@mean_genescore, ($sum_genescore[$j]/$univ)); }

print "\n\tUniv: $univ";
# print "\n\tSum of genescore: @sum_genescore"; #\n\tUniv: $univ";

my @sd_genescore=();
for(my $j=0;$j<=$#colarray;$j++)
{
	my $sqsum=0;
	foreach(@genearray)
	{
		$sqsum+=(${$genescore{$_}}[$j]-$mean_genescore[$j])**2;
	}

	push(@sd_genescore, sqrt($sqsum/$univ));
}


# Indexing gs descriptions and member genes
################################################

$elapsed=0; $elapsed = runtime();
print "\n$elapsed: Indexing gs descriptions and genes ...";

my %gsdesc=();
foreach(@gsdesc)
{
	$_=~s///g;
	@p=(); @p=split("\t",$_);

	$gsdesc{$p[0]}=$p[1];
}

my %gsgenes=(); my %gssize=(); my $gsid;
foreach my $gene (@gset)
{
	$gene=~s///g;
	@p=(); @p=split("\t",$gene);

	if(exists $backgroundlist{$p[0]})
	{
		for(my $j=1;$j<=$#p;$j++)
		{
			$gsid=''; $gsid=$p[$j];
			$gsid=~s/ /_/g;

			push(@{$gsgenes{$gsid}},$p[0]);

			$gssize{$gsid}++;
		}
	}
}

my @gsarray=(); my $totnumgs=0; my $numgswithenoughgenes=0; my $minnumgenes=10;

$totnumgs=scalar(keys %gssize);

foreach my $gs (sort {$gssize{$b}<=>$gssize{$a}} keys %gssize)
{
	if($gssize{$gs}>=$minnumgenes)
	{
		$numgswithenoughgenes++;
		push(@gsarray, $gs);
	}
	else { last; }
}

print "\n\tTotal no. of gss: $totnumgs\n\tNo. of gss with $minnumgenes or more genes: $numgswithenoughgenes\n\n";


# Calculating gene set scores, enrichment z-scores and their pvalues
####################################################################

#print JH "GS.ID\tGS.Desc\tGS.Size";
#@p=(); @p=split("\t",$f[0]);
#for(my $j=0;$j<=$#colarray;$j++)
#{
#	print JH "\t$p[$colarray[$j]]";
#}
#print JH "\n";

$elapsed=0; $elapsed = runtime();
print "$elapsed: Calculating gs scores ...", "\n";

my %gsscore=(); my %gsparzscore=(); my %gsparpval=(); my $numtests=0;

my($meangsscore, $parsd, $zscore, $pvalue, $meantoprint);
for(my $i=0;$i<=$#gsarray;$i++)
{
	$numtests++;
	my $gssize = $gssize{$gsarray[$i]};

	@{$gsscore{$gsarray[$i]}} = @{meangsscore(\@{$gsgenes{$gsarray[$i]}}, $gssize)};	# Mean geneset score

	@{$gsparzscore{$gsarray[$i]}}=();	# Z-score
	@{$gsparpval{$gsarray[$i]}}=();		# P-value

	for(my $j=0;$j<=$#colarray;$j++)
	{
		$meangsscore = ${$gsscore{$gsarray[$i]}}[$j];
		
		$parsd = ($sd_genescore[$j]/sqrt($gssize));

		$zscore = sprintf("%.3f", ($meangsscore - $mean_genescore[$j])/$parsd);
		push(@{$gsparzscore{$gsarray[$i]}}, $zscore);

		$pvalue = sprintf("%.5g", (Statistics::Distributions::uprob (abs($zscore)))*2);
		push(@{$gsparpval{$gsarray[$i]}}, $pvalue);
	
		$meantoprint = sprintf("%.3f", $meangsscore);
	}
}


# Multiple Hypotheses Testing correction
############################################################

$elapsed=0; $elapsed = runtime();
# print "\n$elapsed: Calculating corrected qvalues ...\n";

my %gsesparfdr=(); my($minparpval, $rank, $fdr, $esparfdr);
for(my $j=0;$j<=$#colarray;$j++)
{
	$minparpval=1;
	foreach my $gs (keys %gsparpval)
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
		$esparfdr = (${$gsscore{$_}}[$j]/abs(${$gsscore{$_}}[$j]))*(-1)*log($fdr)/log(10);
		
		push(@{$gsesparfdr{$_}}, $esparfdr);

		$rank--;
	}
}


# Printing results
##################

$elapsed=0; $elapsed = runtime();
print "$elapsed: Printing results ...\n";

my @gsintable = ();	# Array index i: geneset index; j: row in table - gsid, gsdesc, gssize, gsscores
my %gsingraph = ();	# Hash key i: condition; j: geneset pointing to score
foreach my $gs (sort keys %gsesparfdr)
{
	my $inclpar=0;
	for(my $j=0;$j<=$#colarray;$j++)
	{
		if (${$gsesparfdr{$gs}}[$j] >= (-1)*log($qvalcutoff)/log(10))
		{
			$inclpar=1;
			$gsingraph{$colarray[$j]}{$gs}=${$gsparzscore{$gs}}[$j];
		}
	}

	if($inclpar eq 1)
	{
		push @gsintable, [$gs, $gsdesc{$gs}, $gssize{$gs}, @{$gsparzscore{$gs}}];
	}
}

$elapsed=0; $elapsed = runtime();
print "$elapsed: DONE ...\n\n";


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
			$totscore[$k]+=${$genescore{$_}}[$k];
		}
	}

	my @meanscore=();
	for(my $k=0;$k<=$#colarray;$k++) { push(@meanscore, ($totscore[$k]/$size)); }

	return \@meanscore;
}

