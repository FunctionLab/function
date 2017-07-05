$fi1=$ARGV[0];
$fi2=$ARGV[1];
chomp($fi3=$ARGV[2]);

if($#ARGV<2)
{
	print "\nCompares 3 lists and calculates quantities pertaining to different regions of a 3-set Venn diagram\n\n";
	print "<List1> <List2> <list3>\n\n";
	exit;
}

open FH, "$fi1";
chomp(@f=<FH>);
close FH;

open GH, "$fi2";
chomp(@g=<GH>);
close GH;

open HH, "$fi3";
chomp(@h=<HH>);
close HH;

%if=(); %ig=(); %ih=();
foreach(@f){ $_=~s///g; @p=(); @p=split("\t",$_); $if{$p[0]}++; }
foreach(@g){ $_=~s///g; @p=(); @p=split("\t",$_); $ig{$p[0]}++; }
foreach(@h){ $_=~s///g; @p=(); @p=split("\t",$_); $ih{$p[0]}++; }

$s1=scalar(keys %if);$s2=scalar(keys %ig);$s3=scalar(keys %ih);

$common_fg=0; $common_gh=0; $common_fh=0;
$common=0;

foreach(keys %if)
{
	# $v=0;

	if(exists $ig{$_})
	{
		$common_fg++; # $v=1;
		if(exists $ih{$_})
		{
			$common++;
			# print "$_\n";
		}
		# else
		# {
		#	print "$_\n";
		# }
	}

	if(exists $ih{$_})
	{
		$common_fh++; $v=1;
		# unless(exists $ig{$_})
		# {
		#	print "$_\n";
		# }
	}

	# if($v eq 0)
	# {
	#	print "$_\n";
	# }
}

foreach(keys %ig)
{
	if(exists $ih{$_}){ $common_gh++; }
}

$p=$common_fg-$common; $q=$common_fh-$common; $r=$common_gh-$common;
$only_set1 = $s1-($common_fh+$p); $only_set2 = $s2-($common_fg+$r); $only_set3 = $s3-($common_gh+$q);

$union=0; $union=$only_set1+$only_set2+$only_set3+$common_fg+$common_gh+$common_fh+$common;

print "\nNo. of genes in Set1: $s1\nNo. of genes in Set2: $s2\nNo. of genes in Set3: $s3\n\n";
print "No. of genes in the union: $union\nNo. of genes in common: $common\n\n";
print "No. of genes common only to 1 and 2: $p\nNo. of genes common only to 2 and 3: $r\nNo. of genes common only to 1 and 3: $q\n\n";
print "No. of genes in Set1_only: $only_set1\nNo. of genes in Set2_only: $only_set2\nNo. of genes in Set3_only: $only_set3\n\n";
