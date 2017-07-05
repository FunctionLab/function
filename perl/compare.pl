$fi1=$ARGV[0];
chomp($fi2=$ARGV[1]);

open FH, "$fi1";
chomp(@f=<FH>);
close FH;

open GH, "$fi2";
chomp(@g=<GH>);
close GH;

%all=();
%same=();
foreach(@g)
{
	$_=~s///g;
	@p=();
	@p=split("\t",$_);
	
	$same{$p[0]}=$p[1];
	
	$all{$p[0]}++;
}

$s1=scalar(@f);$s2=scalar(@g);
$common=0;

#print "\n";
foreach(@f)
{
	$_=~s///g;
	@p=();
	@p=split("\t",$_);

	$all{$p[0]}++;

	if(exists($same{$p[0]}))
	{
		$common++;
		print "$_\t$same{$p[0]}\n";
	}
#	else{print "$_\n";}
}

$c1=0; $u=0;
foreach(keys %all)
{
	$u++;
	
	if($all{$_}>1)
	{
		$c1++;
	}
}

$only_set1 = $s1-$common; $only_set2 = $s2-$common; $union=$only_set1+$only_set2+$common;
print "\nNo. of genes in Set1: $s1\nNo. of genes in Set2: $s2\nNo. of genes in common: $common\n";
print "No. of genes in Set1_only: $only_set1\nNo. of genes in Set2_only: $only_set2\nTotal no. of genes (Union): $union\n\n";

$c2=0; $c2=$s1+$s2-$u;
#print "Common1: $c1\tCommon2: $c2\nUnion: $u\n\n";

$j=0;
$j=$common/$union;
print "Jaccard = $j\n";

$deno=0; if($s1<$s2){$deno=$s1;} else{$deno=$s2;}
$ovlp=0; $ovlp=$common/$deno;
print "\nMeetMin = $ovlp\n\n";
