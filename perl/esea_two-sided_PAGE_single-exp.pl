#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Time::SoFar qw( runtime );
use Statistics::Distributions;

my ($help, $idat, $igmt, $icover);
my $itail = '2tail'; my $imiss = 'NA'; my $imode = 'both';

pod2usage( -exitstatus => 2, -verbose => 2 ) if ( @ARGV < 2 );
GetOptions( 'help' => \$help,
          'idat=s' => \$idat,
          'igmt=s' => \$igmt,
        'icover=s' => \$icover,
         'itail=s' => \$itail,
         'imiss=s' => \$imiss,
         'imode=s' => \$imode,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

$fout1 = $igmt; $fout1 =~ s/\.\.\///g; $fout1 =~ s/^.*\///g; $fout1 =~s/\.genesets//g; $fout1 =~s/\.gmt//g;
$fout1 = $idat.'-'.$fout1.'.esea.zscore.mat';
$fout2 = $fout1; $fout2 =~ s/\.zscore\.mat$/\.esfdr\.mat/g;
open JH,">$fout1";
open HH,">$fout2";

print "\n\tInput file: $idat\n\tOutput files:\n\t\t$fout1\n\t\t$fout2\n";
if($icover eq 'yes') { $fout3 = $fout1; $fout3 =~ s/\.zscore\.mat$/\.setcover\.mat/g; open SH,">$fout3"; print "\t\t$fout3\n";}


my $time = runtime();
print "\n$time: Indexing input edges, scores and background lists ...\n";

my %escore = (); my @earray = ();
my %bgenes = (); my %bedges = ();
my $mean_escore = 0; my $ntote = 0; my (@p, $e, @tgenes);

open DAT, "$idat" or die "Can't open $idat!";
while (<DAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $e = join '__', sort($p[0], $p[1]);
    push(@earray, $e);
	$bgenes{$p[0]}++; $bgenes{$p[1]}++;
	$bedges{$e}++;

	$escore{$e} = $p[2]; $ntote++;
    $mean_escore += ($p[2] - $mean_escore)/$ntote; }
close DAT;

my $nmisse = 0;
if($imiss eq '0') {
	@tgenes = keys %bgenes;

	for(my $i=0; $i<$#tgenes; $i++) {
		for(my $j=($i+1); $j<=$#tgenes; $j++) {
			$e = join '__', sort($tgenes[$i], $tgenes[$j]);

			unless(exists $bedges{$edge}) {
                $ntote++; $mean_escore += -$mean_escore/$ntote;
				$nmisse++; $bedges{$e}++; } } } }

my $sd_escore = 0;
foreach (@earray) {
	$sd_escore += ( $escore{$_} - $mean_escore )**2; }
if($imiss eq '0') { $sd_escore += ($nmisse*($mean_escore**2)); }

$sd_escore = sqrt($sd_escore/$ntote);

print "\tTotal no. edges: $ntote\n";
print "\tMean escore: $mean_escore\n\tSD escores: $sd_escore\n";

$time = runtime();
print "\n$time: Indexing gs genes, desc and edges ...";

my %gs_genes = (); my %gene_gs = ();
my %gs_gsize = (); my %gs_desc = ();
my ($gs, $desc);

open GMT, "$igmt" or die "Can't open $igmt!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $gs = shift(@p);
    ($desc = shift(@p)) =~ s/ \([0-9]*\)$//g;

    foreach my $g (@p) {
        unless(exists $bgenes{$g}) { next; }
        $gs_genes{$gs}{$g}++;
        $gene_gs{$g}{$gs}++;
        $gs_gsize{$gs}++; } }
close GMT;

my %gs_edges = (); my %gs_esize = ();
my %ann_edges = (); my %coann_edges = ();

foreach my $gs (keys %gs_genes) {
    @tgenes = keys %gs_genes;
	
	for(my $i=0; $i<$#tgenes; $i++) {
		for(my $j=($i+1); $j<=$#tgenes; $j++) {
            $e = join '__', sort($tgenes[$i], $tgenes[$j]);
			
			unless(exists $bedges{$e}) { next; }
            $ann_edges{$e}++;
            $coann_edges{$e}{$gs}++;

            if(($imode eq 'intra') or ($imode eq 'both')) {
                $gs_edges{$gs}{$e}++;
                $gs_esize{$gs}++; } } } }

my $ncoanne = scalar keys %coann_edges;

my $ncranne=0;
foreach my $e (keys %bedges) {
    @p = split '\t', $e;

	if((not exists $coann_edges{$e}) and
        (exists $gene_gs{$p[0]}) and (exists $gene_gs{$p[1]})) {
        $ncranne++;

        if($imode eq 'intra') { next; }
        
        for($k=0;$k<=$#{$gene_gs{$p[0]}};$k++)
        {
            for($l=0;$l<=$#{$gene_gs{$p[1]}};$l++)
            {
                @gsorder = sort(${$gene_gs{$p[0]}}[$k], ${$gene_gs{$p[1]}}[$l]);
                $gs = join '|', @gsorder;
                push(@{$gs_edges{$gs}}, $e);
                $ann_edges{$e}++;

                $gs_esize{$gs}++;
                unless(exists $gspairdesc{$gs})
                {
                    $gspairdesc{$gs}=$gs_desc{$gsorder[0]}.'|'.$gs_desc{$gsorder[1]};
                }
            }
        }
    }
}

@gsarray=(); $totnumgs=0; $numgswithenoughedges=0; $minnumedges=10;

$totnumgs=scalar(keys %gs_esize);
$num_intrags=0; $num_intergs=0;

#%dyn_gs_edges=();
%dyn_set_size=(); %dyn_edge_gs=();
foreach $gs (sort {$gs_esize{$b}<=>$gs_esize{$a}} keys %gs_esize)
{
	if($gs_esize{$gs}>=$minnumedges)
	{
		$numgswithenoughedges++;
		push(@gsarray, $gs);

		$dyn_set_size{$gs}=$gs_esize{$gs};
		foreach $edge (@{$gs_edges{$gs}})
		{
			#$dyn_gs_edges{$gs}{$edge}++;
			$dyn_edge_gs{$edge}{$gs}++;
		}
		
		if($gs=~/\|/)
		{
			$num_intergs++;
			%gs_genes=();
			foreach $edge (@{$gs_edges{$gs}})
			{
				@p=(); @p=split("__",$edge);
				$gs_genes{$p[0]}++; $gs_genes{$p[1]}++;
			}
			$gs_gsize{$gs}=scalar(keys %gs_genes);
		}
		else
		{
			$num_intrags++;
		}
	}
	else { last; }
}

$nanne=scalar(keys %ann_edges);

print "\n\tTotal no. of annotated edges: $nanne\n\tNo. coannotated edges: $ncoanne\n\tNo. cross-annatated edge: $ncranne\n";
print "\n\tTotal no. of gss: $totnumgs\n\tNo. of gss with $minnumedges or more edges: $numgswithenoughedges\n";
print "\tNo. intrags: $num_intrags\n\tNo. intergs: $num_intergs\n\n";

# Calculating edge set scores, enrichment z-scores and their pvalues, and printing z-scores
###########################################################################################

print JH "GS.ID\tDesc\tNo.Genes\tNo.Edges\tSum.Score\tMean.Score\tZ.Score\n";

$time=0; $time = runtime();
print "$time: Calculating gs scores and printing z-scores ...", "\n";

%gsmeanscore=(); %gssumscore=(); %gszscore=();
%gspval=(); $numtests=0;

foreach $gs (@gsarray)
{
	$numtests++;
	$size = $gs_esize{$gs};
	$meanscore = meangsscore(\@{$gs_edges{$gs}}, $size);
	$parsd = ($sd_escore/sqrt($size));
	$zscore = (($meanscore - $mean_escore)/$parsd);

	$pvalue=1;
	if($itail eq 'twotail') { $pvalue = sprintf("%.6g", (Statistics::Distributions::uprob (abs($zscore)))*2); }
	elsif($itail eq 'utail') { $pvalue = sprintf("%.6g", (Statistics::Distributions::uprob ($zscore))); }
	elsif($itail eq 'ltail') { $pvalue = sprintf("%.6g", (1 - Statistics::Distributions::uprob ($zscore))); }
	
	# if($itail eq 'twotail') { $pvalue = (Statistics::Distributions::uprob (abs($zscore))*2); }
	# elsif($itail eq 'utail') { $pvalue = (Statistics::Distributions::uprob ($zscore)); }
	# elsif($itail eq 'ltail') { $pvalue = (1 - Statistics::Distributions::uprob ($zscore)); }
	$gspval{$gs} = $pvalue;

	$gsmeanscore{$gs} = sprintf("%.3f", $meanscore);
	$gssumscore{$gs} = sprintf("%.3f", ($meanscore*$size));
	$gszscore{$gs} = sprintf("%.3f", $zscore);

	print JH "$gs\t$gspairdesc{$gs}\t$gs_gsize{$gs}\t$size";
	print JH "\t$gssumscore{$gs}\t$gsmeanscore{$gs}\t$gszscore{$gs}\n";
}

close JH;

# Multiple Hypotheses Testing correction
########################################

print HH "GS.ID\tDesc\tNo.Genes\tNo.Edges\tSum.Score\tMean.Score\tES.FDR\n";

$time=0; $time = runtime();
print "\n$time: Calculating corrected qvalues ...\n";

%gsesfdr=(); $minpval=1;
foreach $gs (keys %gspval)
{
	if(($gspval{$gs} ne 0)&&($gspval{$gs} < $minpval))
	{
		$minpval=$gspval{$gs};
	}
}

$rank=$numtests;
foreach $gs (sort {$gspval{$b} <=> $gspval{$a}} keys %gspval)
{
	if($gspval{$gs} eq 0) { $gspval{$gs} = 0.1*sprintf("%.6g", $minpval); }

	$fdr = $gspval{$gs}*($numtests/$rank); if($fdr>1) { $fdr = 1; }

	$sign=1; if($gsscore{$gs}<0) { $sign=-1; }
	$gsesfdr{$gs} = $sign*(-1)*log($fdr)/log(10);

	$rank--;
}


# Printing Qvalue-based enrichment scores results
#################################################

$time=0; $time = runtime();
print "\n$time: Printing qvalue-based enrichment scores ...\n";

foreach $gs (keys %gsesfdr)
{
	$esfdrtoprint = sprintf("%.3f", $gsesfdr{$gs});
	print HH "$gs\t$gspairdesc{$gs}\t$gs_gsize{$gs}\t$gs_esize{$gs}";
	print HH "\t$gssumscore{$gs}\t$gsmeanscore{$gs}\t$gszscore{$gs}\t$esfdrtoprint\n";
}

close HH;

if($icover eq 'no')
{
	$time = runtime();
	print "$time: DONE\n\n";
}

else
{
# Set-Cover
###########

	$time = runtime();
	print "\n$time: Performing set-cover and printing results ...\n";

	print SH "GS.ID\tDesc\tNo.Genes\tNo.Edges\tSum.Score\tMean.Score\tZ.Score\tNo.New.Edges\tFrac.Edges.Covered\tCostBenf.Ratio\n";

	%uncovered_edges=%ann_edges; $num_uncovered_edges=scalar(keys %ann_edges);
	%unselected_sets=%gspval; #%selected_sets=();
	%covered_edges=(); $ncovere=0;

	print "\tTot. no. uncovered edges: $num_uncovered_edges\n"; $count=0;
	while($num_uncovered_edges > 0)
	{
		$count++;
		$mincostbenf=1;
		foreach $gs (keys %unselected_sets)
		{
			if($gspval{$gs} eq 1) { $gspval{$gs}=(1-10**(-15));}
			if($dyn_set_size{$gs}>0)
			{
				if((log($gspval{$gs}) eq 0)||($dyn_set_size{$gs} eq 0)) { print "\t\t>> $gs\t$gspval{$gs}\t$dyn_set_size{$gs}\n"; exit; }
				$ratio=(log(10)/((-1)*log($gspval{$gs})*$dyn_set_size{$gs}));
				if($ratio < $mincostbenf) { $mincostbenf=$ratio; $mings=$gs; }
			}
		}

		$ncovere+=$dyn_set_size{$mings};
		$num_uncovered_edges=$nanne-$ncovere;
		
		$mincostbenf_toprint=sprintf("%.3g", $mincostbenf);
		$frac_edges_covered=sprintf("%.3g", ($ncovere/$nanne));

		print SH "$mings\t$gspairdesc{$mings}\t$gs_gsize{$mings}\t$gs_esize{$mings}";
		print SH "\t$gssumscore{$mings}\t$gsmeanscore{$mings}\t$gszscore{$mings}";
		print SH "\t$dyn_set_size{$mings}\t$frac_edges_covered\t$mincostbenf_toprint\n";

		#$selected_sets{$mings}=$mincostbenf;
		delete $unselected_sets{$mings};
	
		foreach $edge (@{$gs_edges{$mings}})
		{
			foreach $gs (keys %{$dyn_edge_gs{$edge}})
			{
				unless($gs eq $mings) { $dyn_set_size{$gs}--; }
			}
		}

		if(($count%500) eq 0)
		{
			$time = runtime();
			print "\t$time:\t$mings\t$gspairdesc{$mings}\t$gs_gsize{$mings}\t$gs_esize{$mings}\n";
			print "\t\tNo. new edges added: $dyn_set_size{$mings}\tCost-Benefit ratio: $mincostbenf_toprint\n";
			print "\t\tFrac. edges covered: $frac_edges_covered\tNo. uncovered edges: $num_uncovered_edges\n";
		}
	}

	$time = runtime();
	print "$time: DONE\n\n";

	close SH;
}


# Subroutine to calculate the mean score of a gs
#####################################################

sub meangsscore
{
	my ($tempga, $size) = @_;
	my $totscore=0; my $meanscore=0;
	
	foreach (@$tempga)
	{
		my $score=0; if(exists $escore{$_}) { $score=$escore{$_}; }
		$totscore+=$score;
	}

	$meanscore=($totscore/$size);
	return $meanscore;
}




__END__

=head1

<Brief_Desc>

=head1 USAGE

./code.pl [--i INPUT_FILE] [--p OPTION] [--o OUTPUT_FILE] [--help]

=head1 DESCRIPTION

This script ...

=head1 ARGUMENTS

=over 12

=item C<--i>

Input file

=item C<--o>

Output file

=item C<--help>

Prints this help message.

=back

=head1 PARAMETERS

=over 12

=item C<--p>

Parameter

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2010 May 15

=cut

