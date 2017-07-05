$f1=$ARGV[0];		# Matrix to normalize; Contains header row.
$nc=$ARGV[1];		# Choice of normalization: "DivColSum", "DivFracNumNonzeroCol", "MulOneminusFracNumNonzeroCol", "DivLogColSum", "DivRowSum", "DivFracNumNonzeroRow", "MulOneminusFracNumNonzeroRow", "DivLogRowSum", "RowCenter", "Binarize:Threshold"
chomp($f3=$ARGV[2]);	# Output: Normalized matrix

open FH,"$f1";
chomp(@f=<FH>);
close FH;
#open GH,"$f2";
#chomp(@g=<GH>);
#close GH;

open HH,">$f3";

print HH "$f[0]\n";
@p=();
@p=split("\t",$f[0]);
$numcol=0; $numcol=$#p;

if($nc=~/^Binarize/) { @p=(); @p=split(":",$nc); $nc=$p[0]; $bint=$p[1];}
@genearray=(); %rowsum=(); %rownumnonzero=(); %colsum=(); %colnumnonzero=();
$totnumgenes=0;
for($i=1;$i<=$#f;$i++)
{
	$totnumgenes++;

	@p=();
	@p=split("\t",$f[$i]);

	push(@genearray, $p[0]);

	if($nc eq 'Binarize') { print HH "$p[0]"; }
	
	for($j=1;$j<=$#p;$j++)
	{
		if($nc eq 'Binarize')
		{
			if($p[$j]>$bint) { print HH "\t1"; }
			else { print HH "\t0"; }
		}
		else
		{
			$mat[$i-1][$j-1]=$p[$j];

			if($nc=~/Row/) { $rowsum{$i-1}+=$p[$j]; if($p[$j] ne 0) { $rownumnonzero{$i-1}++; } }
			if($nc=~/Col/) { $colsum{$j-1}+=$p[$j]; if($p[$j] ne 0) { $colnumnonzero{$j-1}++; } }
		}
	}

	if($nc eq 'Binarize') { print HH "\n"; }
}

unless($nc eq 'Binarize')
{
	for($i=0;$i<=$#genearray;$i++)
	{
		print HH "$genearray[$i]";

		$rowmean=0; $rowmean=($rowsum{$i}/$numcol);
		for($j=0;$j<$numcol;$j++)
		{
			$val=0;

			if($nc eq 'DivFracNumNonzeroCol') { $val=($mat[$i][$j]/$colnumnonzero{$j})*$totnumgenes; }
			elsif($nc eq 'MulOneminusFracNumNonzeroCol') { $val=$mat[$i][$j]*(1-($colnumnonzero{$j}/$totnumgenes)); }
			elsif($nc eq 'DivColSum') { $val=($mat[$i][$j]/$colsum{$j}); }
			elsif($nc eq 'DivLogColSum') { $val=($mat[$i][$j]/(log($colsum{$j})/log(2))); }

			elsif($nc eq 'DivFracNumNonzeroRow') { $val=($mat[$i][$j]/$rownumnonzero{$i})*$numcol; }
			elsif($nc eq 'MulOneminusFracNumNonzeroRow') { $val=$mat[$i][$j]*(1-($rownumnonzero{$i}/$numcol)); }
			elsif($nc eq 'DivRowSum') { $val=($mat[$i][$j]/$rowsum{$i}); }
			elsif($nc eq 'DivLogRowSum') { $val=($mat[$i][$j]/(log($rowsum{$i})/log(2))); }

			elsif($nc eq 'RowCenter') { $val=$mat[$i][$j]-$rowmean; }

			$printval=0; $printval=sprintf("%.3f", $val);
			print HH "\t$printval";
		}

		print HH "\n";
	}
}

close HH;
