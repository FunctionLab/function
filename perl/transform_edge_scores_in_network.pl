$f1=$ARGV[0];		# Network file: DAT format
$t=$ARGV[1];		# Transformation parameter: '0to1', '-1to1', 'inv-Fisher'
chomp($f3=$ARGV[2]);	# Outfile: Network in DAT format

open FH,"$f1";
chomp(@f=<FH>);
close FH;

open HH,">$f3";

if($t eq '0to1') { $maxscore=0; }
foreach(@f)
{
	$_=~s///g;

	@p=();
	@p=split("\t",$_);

	if($t eq '0to1') { if($p[2]>$maxscore) { $maxscore=$p[2]; } }
	else
	{
		$tscore=0;
		if($t eq '-1to1') { $tscore=((2*$p[2])-1); }
		elsif($t eq 'inv-Fisher') { $tscore=sprintf("%.6f",((exp($p[2])-1)/(exp($p[2])+1))); }

		print HH "$p[0]\t$p[1]\t$tscore\n";
	}
}

if($t eq '0to1')
{
	foreach(@f)
	{
		@p=();
		@p=split("\t",$_);

		$tscore=($p[2]/$maxscore);
		print HH "$p[0]\t$p[1]\t$tscore\n";
	}
}

close HH;
