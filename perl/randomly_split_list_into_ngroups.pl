#!/usr/bin/perl
use strict;
use warnings;
use List::Util 'shuffle';

$f1=$ARGV[0];		# List of entities; Annotations are allowed in subsequent columns
$ngrp=$ARGV[1];		# No. of groups to split into
chomp($f3=$ARGV[2]);	# Random group indexed list

open FH,"$f1";
chomp(@f=<FH>);
close FH;
#open GH,"$f2";
#chomp(@g=<GH>);
#close GH;

open HH,">$f3";

@origarray=(); %ann=();
foreach(@f)
{
	$_=~s///g;

	@p=();
	@p=split("\t",$_);

	push(@origarray, $p[0]);

	if($#p>0)
	{
		for($j=1;$j<=$#p;$j++)
		{
			push(@{$ann{$p[0]}}, $p[$j]);
		}
	}
}

$totnum=0; $totnum=scalar(@origarray);
$numineachgrp=int(($totnum/$ngrp)+0.5);

print "\nTotal no. entities: $totnum\nNo. in each of the $ngrp groups: $numineachgrp\n";

@newarray=();
@newarray=shuffle(@origarray);
#print "\n@origarray\n@newarray\n\n";

$grpidx=1;
for($i=0;$i<=$#newarray;$i++)
{
	if($grpidx eq $ngrp)
	{
		do
		{
			print HH "$grpidx\t$newarray[$i]";
			
			if(exists $ann{$newarray[$i]})
			{
				for($j=0;$j<=$#{$ann{$newarray[$i]}};$j++)
				{
					print HH "\t${$ann{$newarray[$i]}}[$j]";
				}
			}

			print HH "\n";
			
			$i++;
		}while($i<=$#newarray);
	}
	else
	{
		$numinthisgrp=0;
		do
		{
			$numinthisgrp++;

			print HH "$grpidx\t$newarray[$i]";
			
			if(exists $ann{$newarray[$i]})
			{
				for($j=0;$j<=$#{$ann{$newarray[$i]}};$j++)
				{
					print HH "\t${$ann{$newarray[$i]}}[$j]";
				}
			}

			print HH "\n";

			$i++;
		}while($numinthisgrp<$numineachgrp);

		$grpidx++;
	}

	$i--;
}

close HH;

