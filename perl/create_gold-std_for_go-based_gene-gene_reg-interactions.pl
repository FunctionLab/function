#!/usr/bin/perl
use strict;
use warnings;
# use lib '/Genomics/Users/arjunk/src/perl/lib';
use Getopt::Long;
use Pod::Usage;
use PDF;
use Data::Dumper;
#use Time::SoFar qw(runtime);

# Reading in files and options
# ############################
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);

my($help, $gafile, $regtarpairsfile);
my($out_goldstd, $out_regtarpairs_ann, $out_reggs, $out_coregs);
my $mingenes = 5; my $maxgenes = 100; my $minedges = 30; my $pvalovlp = 0.05;

GetOptions(	'help' => \$help,
	       'iga=s' => \$gafile,
	   'iregtar=s' => \$regtarpairsfile,
	  'pvalovlp=f' => \$pvalovlp,
	     'ogstd=s' => \$out_goldstd,
	   'oregtar=s' => \$out_regtarpairs_ann,
        'oreggs=s' => \$out_reggs,
       'ocoregs=s' => \$out_coregs,
	  'mingenes=i' => \$mingenes,
	  'maxgenes=i' => \$maxgenes,
	  'minedges=i' => \$minedges) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open GA, "$gafile" or die "Can't open file: $gafile!"; chomp(my @ga=<GA>); close GA;
open RT, "$regtarpairsfile" or die "Can't open file: $regtarpairsfile!"; chomp(my @rt=<RT>); close RT;
open GS, ">$out_goldstd"; open ART, ">$out_regtarpairs_ann";
open RGS, ">$out_reggs"; open CGS, ">$out_coregs";

# Associating genes to GOterms and vice-versa
# ###########################################
my @line;
my %gs_size = (); my %gs_gene = (); my %gene_gs = ();
foreach(@ga) {
	@line=split("\t",$_);
	$gs_gene{$line[0]}{$line[1]}++; $gs_size{$line[0]}++;
	$gene_gs{$line[1]}{$line[0]}++;
}

# TryOut: For every reg-tar goterm pair, removing a gene from tar if it is
# annotated to reg
my ($r, $t);
foreach(@rt) {
    @line = split('\t', $_);
    $r = $line[1]; $t = $line[3];

    #if(($line[0]=~/R-C/) and (exists $gs_incl{$r}) and (exists $gs_incl{$t})) {
    if($line[0]=~/R-C/) {
        foreach my $g1 (keys %{$gs_gene{$r}}) {
            if(exists $gs_gene{$t}{$g1}) { delete $gs_gene{$t}{$g1}; }
        }
    }
}

# Indexing coannotated gene pairs
# ###############################
my %gs_incl = (); my %coann_genepair = ();
foreach my $gs (keys %gs_size) {
	if(($gs_size{$gs} >= $mingenes)&&($gs_size{$gs} <= $maxgenes)) {
		$gs_incl{$gs}++;
		my @temp_gs = (keys %{$gs_gene{$gs}});
		for(my $i=0;$i<$#temp_gs;$i++) {
			for(my $j=($i+1);$j<=$#temp_gs;$j++) {
				my $edge = join '__', sort($temp_gs[$i], $temp_gs[$j]);
				$coann_genepair{$edge}++;
			}
		}
	}
}

# Populating regulatory-target GOterm pairs and gene pairs
# ########################################################
print ART "Reg.GOID\tNo.Genes\tGO.Term\tTar.GOID\tNo.Genes\tGO.Term\tNo.Common.Genes\tNo.Gene.Pairs\n";

my ($sr, $st);
my %gs_term =(); my %reg_gs =(); my %core_gs = ();
my %regtar_genepair = (); my %regtar_gspair = ();
my $tot_regtargspairs = 0; my $num_regtargspairs_withenoughgenepairs = 0;
my %gstd_regtar_genepair = ();

foreach(@rt) {
	@line=split("\t",$_);
	$r = $line[1]; $t = $line[3];

    if(($line[0]=~/R-C/) and (exists $gs_incl{$r}) and (exists $gs_incl{$t})) {
        $tot_regtargspairs++;

        $sr = $gs_size{$r}; $st = $gs_size{$t};
        $gs_term{$r} = $line[2]; $gs_term{$t} = $line[4];

        $reg_gs{$r}{$t}++; $core_gs{$t}{$r}++;

        my $num_genepairs = 0;
        foreach my $g1 (keys %{$gs_gene{$r}}) {
            foreach my $g2 (keys %{$gs_gene{$t}}) {
                my $edge = join '__', sort($g1, $g2);
                unless (exists $coann_genepair{$edge}) {
                    $regtar_genepair{$edge}++; $gstd_regtar_genepair{$edge}=1;
                    $num_genepairs++;
                }
            }
        }

        if($num_genepairs >= $minedges) {
            $num_regtargspairs_withenoughgenepairs++;
            $regtar_gspair{$r.'_'.$t}=$line[0];
            print ART "$r\t$sr\t$line[2]\t$t\t$st\t$line[4]\t$num_genepairs\n";
        }
    }
}

# Counting regulatory and core genes
# ##################################
my ($num_reggs, $num_coregs) = (0, 0);
my %reg_genes = (); my %core_genes = ();

foreach my $r (keys %reg_gs) {
    $num_reggs++;
    print RGS "$r\t$gs_term{$r}";

    foreach my $g (keys %{$gs_gene{$r}}) {
        $reg_genes{$g}++;
        print RGS "\t$g";
    }
    print RGS "\n";
}

foreach my $c (keys %core_gs) {
    $num_coregs++;
    print CGS "$c\t$gs_term{$c}";

    foreach my $g (keys %{$gs_gene{$c}}) {
        $core_genes{$g}++;
        print CGS "\t$g";
    }
    print CGS "\n";
}

my $num_reggenes = scalar(keys %reg_genes); my $num_coregenes = scalar(keys %core_genes);

# Comparing GOterms (reg with reg, and, core with core) to calculate overlapping sets
# ###################################################################################
my @reg_gsarray = sort(keys %reg_gs); my %ovlp_reggs = ();
for(my $i=0;$i<$#reg_gsarray;$i++) {
	for(my $j=($i+1);$j<=$#reg_gsarray;$j++) {
		my $common = 0;
		foreach my $g1 (keys %{$gs_gene{$reg_gsarray[$i]}}) {
			foreach my $g2 (keys %{$gs_gene{$reg_gsarray[$j]}}) {
				if($g1 eq $g2) { $common++; }
			}
		}
		my $prob = hypergeometric_tail(scalar(keys %reg_genes), $gs_size{$reg_gsarray[$i]}, $gs_size{$reg_gsarray[$j]}, $common);
		if($prob <= $pvalovlp) { $ovlp_reggs{join('__', $reg_gsarray[$i], $reg_gsarray[$j])}++; }
	}
}

my @core_gsarray = sort(keys %core_gs); my %ovlp_coregs = ();
for(my $i=0;$i<$#core_gsarray;$i++) {
	for(my $j=($i+1);$j<=$#core_gsarray;$j++) {
		my $common = 0;
		foreach my $g1 (keys %{$gs_gene{$core_gsarray[$i]}}) {
			foreach my $g2 (keys %{$gs_gene{$core_gsarray[$j]}}) {
				if($g1 eq $g2) { $common++; }
			}
		}
		my $prob = hypergeometric_tail(scalar(keys %core_genes), $gs_size{$core_gsarray[$i]}, $gs_size{$core_gsarray[$j]}, $common);
		if($prob <= $pvalovlp) { $ovlp_coregs{join '__', $core_gsarray[$i], $core_gsarray[$j]}++; }
	}
}

# Getting negative reg-core gene pairs
# ####################################
my %neg_regtar_genepair = (); my $ovlppar;
foreach my $r (keys %reg_gs) {
	foreach my $t (keys %core_gs) {
		unless(exists $regtar_gspair{$r.'_'.$t}) {
			$ovlppar=0;
			foreach my $rtar (keys %{$reg_gs{$r}}) {
				if(exists $ovlp_reggs{join '__', sort($t, $rtar)}) { $ovlppar=1; last; }
			}
			
			if($ovlppar==0) {
				foreach my $treg (keys %{$core_gs{$t}}) {
					if(exists $ovlp_reggs{join '__', sort($r, $treg)}) { $ovlppar=1; last; }
				}
			}

			if($ovlppar==0) {
				foreach my $g1 (keys %{$gs_gene{$r}}) {
					foreach my $g2 (keys %{$gs_gene{$t}}) {
						my $edge = join '__', sort($g1, $g2);
						unless(($g1 eq $g2)||(exists $coann_genepair{$edge})||(exists $regtar_genepair{$edge})) {
							$neg_regtar_genepair{$edge}++; $gstd_regtar_genepair{$edge}=0;
						}
					}
				}
			}
		}
	}
}

# Recording complete gold-std and printing it
# ###########################################
my @genes; my $num_posgenepairs = 0; my $num_neggenepairs = 0;
foreach my $edge (keys %gstd_regtar_genepair) {
	@genes = split('__', $edge);
	print GS "$genes[0]\t$genes[1]\t$gstd_regtar_genepair{$edge}\n";
	if($gstd_regtar_genepair{$edge} == 1) { $num_posgenepairs++; }
	else { $num_neggenepairs++; }
}

close GS; close ART;
close RGS; close CGS;

print "\nTot. no. reg-tar gs pairs: $tot_regtargspairs\nNo. reg-tar pairs with enough edges: $num_regtargspairs_withenoughgenepairs\n";
print "\nNo. regulatory gs: $num_reggs\nNo. regulatory genes: $num_reggenes\n";
print "\nNo. core gs: $num_coregs\nNo. core genes: $num_coregenes\n";
print "\nNo. positive genepairs: $num_posgenepairs\nNo. negative genepairs: $num_neggenepairs\n\n";

__END__

=head1 NAME

code.pl
	- <Brief_Desc>

=head1 USAGE

./code.pl [--i INPUT_FILE] [--p OPTION] [--o OUTPUT_FILE] [--help]

=head1 DESCRIPTION

This script ...

=head1 ARGUMENTS

=over 12

=item C<--i>

Input file

=item C<--p>

Option

=item C<--o>

Output file

=item C<--help>

prints this documentation

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

=cut
