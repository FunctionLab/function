#################################################################
# The following code converts a matrix of gene expression val-  #
# ues for all the probes in a microarray (obtained from R after #
# normalization) into expression data in PCL format where the   #
# probeids have been mapping to gene ids. This is done assuming #
# that there is one-to-one mapping between the probes and the   #
# genes. If not, use the 'map_matrix_from_one_idspace_to_anoth- #
# er.pl script that does not assume anyhting about the mapping. #
#################################################################


$f1=$ARGV[0];		# Matrix space/tab seperated
$f2=$ARGV[1];		# ID mapping file: <id1> <id2>
chomp($f3=$ARGV[2]);	# PCL format file

open FH,"$f1";
chomp(@f=<FH>);
close FH;
open GH,"$f2";
chomp(@g=<GH>);
close GH;

open HH,">$f3";

%idmap=();
foreach(@g)
{
	$_=~s///g;

	@p=();
	@p=split("\t",$_);

	$idmap{$p[0]}=$p[1];
}

$f[0]=~s///g;
#$f[0]=~s/^ /Probe /g;
@p=();
@p=split(" ",$f[0]);

$numcolumns=0;
print HH "Probe\tGene\tGWEIGHT";
for($j=1;$j<=$#p;$j++)
{
	print HH "\t$p[$j]";
	$numcolumns++;
}
print HH "\nEWEIGHT\t\t";
for($j=1;$j<=$#p;$j++)
{
	print HH "\t1";
}
print HH "\n";

#%usedid2=();
$numgenes=0; # %numcol=();
for($i=1;$i<=$#f;$i++)
{
	$f[$i]=~s///g;

	@p=();
	@p=split(" ",$f[$i]);

	if(exists $idmap{$p[0]})
	{
		# unless(exists $usedid2{$idmap{$p[0]}})
		{
			print HH "$p[0]\t$idmap{$p[0]}\t1";

			for($j=1;$j<=$#p;$j++)
			{
				print HH "\t$p[$j]";
			}

			print HH "\n";

			# $usedid2{$idmap{$p[0]}}++;
			$numgenes++;
			# $numcol{$#p+1}++;
		}
	}
}

print "\n> $f1\nNo. of columns: $numcolumns","; genes: $numgenes\n\n";

#foreach(keys %numcol)
#{
#	print "$_\n";
#}
#print "\n";

close HH;
