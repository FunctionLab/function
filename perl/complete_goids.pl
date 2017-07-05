$f1=$ARGV[0];		# Genes and GOIDs without the GO:0*
chomp($f3=$ARGV[1]);	# Completed GOIDs

open FH,"$f1";
chomp(@f=<FH>);
close FH;

open HH,">$f3";

foreach(@f)
{
	@p=();
	@p=split("\t",$_);

	$goid=$p[1];

	while(length($goid) < 7) { $goid='0'.$goid; }

	$goid='GO:'.$goid;

	print HH "$p[0]\t$goid\n";
}

close HH;
