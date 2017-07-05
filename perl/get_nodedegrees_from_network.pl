$f1=$ARGV[0];		# Network in DAT format
$t=$ARGV[1];		# Type of score transformation: 'abs', 'log10'
chomp($f3=$ARGV[2]);	# Outfile: <Genes> <NodeDegree>

open FH,"$f1";
chomp(@f=<FH>);
close FH;

open HH,">$f3";

%node_degree=();
foreach(@f)
{
	@p=();
	@p=split("\t",$_);

	$tscore=$p[2];
	if($t eq 'abs') { $tscore=abs($p[2]); }
	elsif($t eq 'log10') { $tscore=log($p[2])/log(10); }

	$node_degree{$p[0]}+=$tscore;
	$node_degree{$p[1]}+=$tscore;
}

foreach(sort {$node_degree{$b} <=> $node_degree{$a}} keys %node_degree)
{
	print HH "$_\t$node_degree{$_}\n";
}

close HH;
