#####################################################################################
# This script calculates the level of enrichment of a bunch genesets in the given   #
# list of genes. The given list can be a single list of genes, or can include a 2nd #
# column of 'factor' information, e.g. 'Up', 'Down', 'ClusterID', 'cisElementID',   #
# or 'PathwayID'. If factors are included, enrichment of the genesets will be calc- #
# ulated for each factored-genelist. Enrichment is calculated using the cummulative #
# Hypergeometric test, and (-1)*log(FDR) values are reported as scores.             #
#                                                                                   #
# All the genesets to use are to be with filenames ending in .genesets.             #
# The background list of genes should be the total list of all the genes with anno- #
# tations in the geneset under consideration.                                       #
#                                                                                   #
# Consider this script as a generic tool to compare a pairs of geneset datasets.    #
#####################################################################################

use PDF;
use Time::SoFar qw( runtime );

$max_fdr = 0.01; $min_num_common_genes = 3;
$min_esfdr=0; $min_esfdr=(-1)*(log($max_fdr)/log(10));

my($f1, $f2, $f3) = @ARGV;
#$f1: List of given (test) genes
#$f2: Geneset data: <Gene> <GeneSet>
#$f3: Outfile: Factor-GeneSet pairs with significant overlap

open FH, "$f1" or die "Can't open $f1!"; chomp(my @fac=<FH>); close FH;
open GH, "$f2" or die "Can't open $f2!"; chomp(my @gset=<GH>); close GH;
open JH, ">$f3";

$time = runtime();
print "$time: Indexing test gene and genesets desc and background lists ...";

%fac_desc=(); %fac_bg=(); $tot_facbg=0;
if($f1=~ /.gmt$/) {
    foreach (@fac) {
        $_=~ s///g; @p=(); @p=split("\t",$_);
        $fac_desc{$p[0]} = $p[1];
        for($j=2; $j<=$#p; $j++) {
            $fac_bg{$p[$j]}++;
        }
    }
    $tot_facbg = scalar(keys %fac_bg);
}
else {
    $fac_desc_file=''; $fac_desc_file=$f1.'.desc';
    $fac_bgfile=''; $fac_bgfile=$f1.'.bg';
    open DH1,"$fac_desc_file"; chomp(@facdesc=<DH1>); close DH1;
    open BH1,"$fac_bgfile"; chomp(@facbg=<BH1>); close BH1;

    foreach(@facdesc) {
        $_=~s///g; @p=(); @p=split("\t",$_);
        $fac_desc{$p[0]}=$p[1];
    }

    foreach(@facbg) { $_=~s///g; $tot_facbg++; $fac_bg{$_}++; }
}

%gs_desc=(); %gs_bg=(); $tot_gsbg=0; %bg_list=(); $univ=0;
if($f2 =~ /.gmt$/) {
    foreach (@gset) {
        $_=~ s///g; @p=(); @p=split("\t",$_);
        $gs_desc{$p[0]} = $p[1];
        for($j=2; $j<=$#p; $j++) {
            $gs_bg{$p[$j]}++;
            if(exists $fac_bg{$p[$j]}) { $bg_list{$p[$j]}++; $univ++; }
        }
    }
    $tot_gsbg = scalar(keys %gs_bg);
}
else {
    $gs_desc_file=''; $gs_desc_file=$f2.'.desc';
    $gs_bgfile=''; $gs_bgfile=$f2.'.bg';
    open DH2,"$gs_desc_file"; chomp(@gsdesc=<DH2>); close DH2;
    open BH2,"$gs_bgfile"; chomp(@gsbg=<BH2>); close BH2;

    foreach(@gsdesc) {
        $_=~s///g; @p=(); @p=split("\t",$_);
        $gs_desc{$p[0]}=$p[1];
    }

    foreach(@gsbg) {
        $_=~s///g; $tot_gsbg++; $gs_bg{$_}++;
        if(exists $fac_bg{$_}) { $bg_list{$_}++; $univ++; }
    }
}

print "\n\tTot. no. TG genes: $tot_facbg\n\tTot. no. GS genes: $tot_gsbg\n\tNo. of genes in background list: $univ\n";

%fac_genes=(); %fac_size=(); $tot_fac_genes=0; $num_fac_genes_inbg=0;
$num_fac=0; @fac_array=();
if($f1 =~ /.gmt$/) {
    %temp_genes = (); %temp_genes_inbg;
    foreach (@fac) {
        $_=~s///g; @p=(); @p=split("\t",$_);
        $num_fac++; push(@fac_array, $p[0]);

        for($j=2; $j<=$#p; $j++) {
            $temp_genes{$p[$j]}++;
            if(exists $bg_list{$p[$j]}) {
                $temp_genes_inbg{$p[$j]}++;
                push(@{$fac_genes{$p[0]}}, $p[$j]);
                $fac_size{$p[0]}++;
            }
        }
    }
    
    $tot_fac_genes = scalar(keys %temp_genes);
    $num_fac_genes_inbg = scalar(keys %temp_genes_inbg);
}
else {
    foreach(@fac) {
        $tot_fac_genes++;

        $_=~s///g;
        @p=(); @p=split("\t",$_);

        if(exists $bg_list{$p[0]}) {
            $num_fac_genes_inbg++; $fac='';

            if($#p>0) {
                for($j=1; $j<=$#p; $j++) {
                    $fac=$p[$j];
                    push(@{$fac_genes{$fac}}, $p[0]);
                    $fac_size{$fac}++;
                }
            }
            else {
                $fac=$f1;
                push(@{$fac_genes{$fac}}, $p[0]);
                $fac_size{$fac}++;
            }


        }
    }

    foreach $fac (keys %fac_genes) {
        $num_fac++;
        push(@fac_array, $fac);
    }
}


print "\n\tTotal no. TGs: $tot_fac_genes\n\tNo. TGs in background: $num_fac_genes_inbg\n\tNo. TG factors: $num_fac\n";

%gs_genes=(); %gs_size=();
if($f2 =~ /.gmt$/) {
    foreach (@gset) {
        $_=~s///g; @p=(); @p=split("\t",$_);
        $p[0] =~ s/ /_/g;

        %temp_genes = (); %temp_genes_inbg = ();
        for($j=2; $j<=$#p; $j++) {
            $temp_genes{$p[$j]}++;
            if(exists $bg_list{$p[$j]}) {
                $temp_genes_inbg{$p[$j]}++;
                push(@{$gs_genes{$p[0]}}, $p[$j]);
                $gs_size{$p[0]}++;
            }
        }
    }
}
else {
    foreach (@gset) {
        $_=~s///g; @p=(); @p=split("\t",$_);

        if(exists $bg_list{$p[0]}) {
            for($j=1; $j<=$#p; $j++) {
                $p[$j]=~s/ /_/g;
                push(@{$gs_genes{$p[$j]}},$p[0]);
                $gs_size{$p[$j]}++;
            }
        }
    }
}

$tot_gs = scalar(keys %gs_genes);
@gs_array=(); $num_gs_wenoughgenes=0;
foreach $gs (sort {$gs_size{$b}<=>$gs_size{$a}} keys %gs_size) {
    if($gs_size{$gs}>=3) {
        $num_gs_wenoughgenes++;
        push(@gs_array, $gs);
    }
    else { last; }
}

print "\n\tTotal no. of genesets: $tot_gs\n\tNo. of genesets with 3 or more genes: $num_gs_wenoughgenes\n\n";

$time = runtime();
print "$time: Comparing to genesets to testgenes ...";

$num_tests=0;
%facgs_pval=(); %facgs_num_common=(); # %facgsuniongenenum=();
for($i=0; $i<=$#fac_array; $i++) {
	%temp_fac_genes=(); $sizea=0;
	foreach $facgene (@{$fac_genes{$fac_array[$i]}}) {
		$temp_fac_genes{$facgene}++;
		$sizea++;
	}

	for($j=0; $j<=$#gs_array; $j++) {
        if($fac_array[$i] eq $gs_array[$j]) { next; }
        
        $facgs_tag=''; $facgs_tag=$fac_array[$i].'__'.$gs_array[$j];
        $revfacgs_tag=''; $revfacgs_tag=$gs_array[$j].'__'.$fac_array[$i];
        if(exists $facgs_num_common{$revfacgs_tag}) { next; }
        
        $sizeb=0; $common=0;
        foreach $gsgene (@{$gs_genes{$gs_array[$j]}}) {
            $sizeb++;
            if(exists $temp_fac_genes{$gsgene}) {
                $common++;
            }
        }

        $facgs_num_common{$facgs_tag}=$common;

        if($common>=$min_num_common_genes) {
            $num_tests++;

            $prob=1;
            $prob=hypergeometric_tail($univ, $sizea, $sizeb, $common);
            if($prob le 0) { $prob=1E-200; } elsif($prob>1) { $prob=1; }
            $facgs_pval{$facgs_tag}=$prob;
        }
    }
}

%facgs_fdr=(); $rank=$num_tests;
foreach (sort {$facgs_pval{$b} <=> $facgs_pval{$a}} keys %facgs_pval) {
	$facgs_fdr{$_}=($facgs_pval{$_}*$num_tests)/$rank;
	if($facgs_fdr{$_}>1){ $facgs_fdr{$_}=1; }

	$rank--;
}

$time = runtime();
print "\n\n$time: Printing significant pairs of factorgenes-gs_genes ...";

print JH "Factor\tFac.No.Genes\tFac.Desc\t";
print JH "GeneSet\tGS.No.Genes\tGS.Desc\t";
print JH "No.Common.Genes\tOverlap\tJaccard\tLOR\tES.qval\n";

%fac_wsiggs=(); %gs_wsigfac=();
for($i=0; $i<=$#fac_array; $i++) {
	$sizea=0; $sizea=$fac_size{$fac_array[$i]};

	for($j=0; $j<=$#gs_array; $j++) {
		$facgs_tag=''; $facgs_tag=$fac_array[$i].'__'.$gs_array[$j];

		$sizeb=0; $sizeb=$gs_size{$gs_array[$j]};
		
		$common=0; $common=$facgs_num_common{$facgs_tag};
		
		$min=$sizea; if($sizeb<$min) { $min=$sizeb; }
		$ovlp=0; $ovlp=($common/$min);

		$union=0; $union=$sizea+$sizeb-$common;
		$jac=0; $jac=($common/$union);

        $lor=0;
        unless($common == 0) { $lor = log(($common/$sizeb)/($sizea/$univ))/log(2); }

        #$prec=0; $prec=($common/$sizea);
        #$recall=0; $recall=($common/$sizeb);
        #$pdotr=0; $pdotr=($prec*$recall);
        #$f1_score=0; $f1_score=(2*$pdotr)/($prec+$recall);
        #$fpt5_score=0; $fpt5_score=(5*$pdotr)/($prec+(4*$recall));

		if(exists $facgs_pval{$facgs_tag}) {
			$espval=0; $espval=(-1)*log($facgs_pval{$facgs_tag})/log(10);
			$esqval=0; $esqval=(-1)*log($facgs_fdr{$facgs_tag})/log(10);

			$facgs_esqval{$facgs_tag}=$esqval;
			$facgs_ovlp{$facgs_tag}=$ovlp;

			if($esqval>=$min_esfdr) {
				$fac_wsiggs{$fac_array[$i]}++; $gs_wsigfac{$gs_array[$j]}++;
				
				print JH "$fac_array[$i]\t$sizea\t$fac_desc{$fac_array[$i]}\t";
				print JH "$gs_array[$j]\t$sizeb\t$gs_desc{$gs_array[$j]}\t";
				print JH "$common\t$ovlp\t$jac\t$lor\t$esqval\n"; #\t$pdotr\t$f1_score\t$fpt5_score\n";
			}
		}
	}
}

$num_fac_wsiggs=0;
foreach $fac (keys %fac_wsiggs) { $num_fac_wsiggs++; }

$num_gs_wsigfac=0;
foreach $gs (keys %gs_wsigfac) { $num_gs_wsigfac++; }

print "\nNo. factors with significant geneset enrichment: $num_fac_wsiggs\n";
print "No. genesets with significant factor enrichment: $num_gs_wsigfac\n";

close JH;
$time = runtime();
print "$time:: DONE\n\n";

