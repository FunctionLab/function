$f1=$ARGV[0];		# List of indexed genes: <Index> <Gene> ...
#$f2=$ARGV[1];
chomp($net=$ARGV[1]);	# Gold-standard network to be split for CV

open FH,"$f1";
chomp(@f=<FH>);
close FH;
#open GH,"$f2";
#chomp(@g=<GH>);
#close GH;

#open HH,">$f3";
#close HH;

%group_genes=(); %all_genes=();
foreach(@f)
{
	@p=();
	@p=split("\t",$_);

	$all_genes{$p[1]}++;
	push(@{$group_genes{$p[0]}}, $p[1]);
}

$numfold=scalar(keys %group_genes);
print "$numfold folds ...\n";
# $numfold.'-fold_cv.'.'train.fold12'
# 3-fold_cv.train.12.genes; 3-fold_cv.test.3.genes

$tag_fold = join("", keys %group_genes);
$tag_train = $numfold.'-fold_cv.train.';
$tag_test = $numfold.'-fold_cv.test.';

foreach $f (sort {$a <=> $b} keys %group_genes)
{
	$tag_train_fold=$tag_fold; $tag_train_fold=~s/$f//g;
	$tag_test_fold=$f;

	%temp_genes=();
	
	$teo=$tag_test.$tag_test_fold.'.genes';
	open TEG, ">$teo";
	foreach $g (@{$group_genes{$f}})
	{
		$temp_genes{$g}++;		
		print TEG "$g\n";
	}
	close TEG;

	$tro=$tag_train.$tag_train_fold.'.genes';
	open TRG, ">$tro";
	foreach $g (keys %all_genes)
	{
		unless(exists $temp_genes{$g})
		{
			print TRG "$g\n";
		}
	}
	close TRG;
}
