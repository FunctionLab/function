#!/usr/bin/perl
use strict;
use warnings;
# use lib '/Genomics/Users/arjunk/src/perl/lib';
use Getopt::Long;
use Pod::Usage;
use PDF;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);

my($help, $gafile, $regtarpairsfile, $out_goldstd, $out_regtarpairs_ann);
my $mingenes = 5; my $maxgenes = 100; my $minedges = 10; my $ovlppval = 0.05;
GetOptions(	  'help' => \$help,
		   'i=s' => \$gafile,
		'ireg=s' => \$regtarpairsfile,
		'ovlp=f' => \$ovlppval,
		 'ogs=s' => \$out_goldstd,
		'oreg=s' => \$out_regtarpairs_ann,
	    'mingenes=i' => \$mingenes,
	    'maxgenes=i' => \$maxgenes,
	    'minedges=i' => \$minedges) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open GA, "$gafile" or die "Can't open file: $gafile!"; chomp(my @ga=<GA>); close GA;
open RT, "$regtarpairsfile" or die "Can't open file: $regtarpairsfile!"; chomp(my @rt=<RT>); close RT;
open GS, ">$out_goldstd"; open ART, ">$out_regtarpairs_ann";

my @line;
my %gs_size = (); my %gs_gene = (); my %gene_gs = ();
foreach(@ga) {
	@line=split("\t",$_);
	$gs_gene{$line[0]}{$line[1]}++; $gs_size{$line[0]}++;
	$gene_gs{$line[1]}{$line[0]}++;
}

my %gs_considered = (); my %coann_genepair = ();
foreach my $gs (keys %gs_size) {
	if(($gs_size{$gs} >= $mingenes)&&($gs_size{$gs} <= $maxgenes)) {
		$gs_considered{$gs}++;
		my @temp_gs = (keys %{$gs_gene{$gs}});
		for(my $i=0;$i<$#temp_gs;$i++) {
			for(my $j=($i+1);$j<=$#temp_gs;$j++) {
				my $edge = join '__', sort($temp_gs[$i], $temp_gs[$j]);
				$coann_genepair{$edge}++;
			}
		}
	}
}

my ($r, $t, $sr, $st);
my %gs_term =(); my %reg_gs =(); my %core_gs = ();
my %regtar_genepair = (); my %regtar_gspair = ();
my $tot_regtargspairs = 0; my $num_regtargspairs_withenoughgenepairs = 0;
foreach(@rt) {
	@line=split("\t",$_);
	$r = $line[1]; $t = $line[3];

	if((exists $gs_considered{$r})&&(exists $gs_considered{$t})) {
		$tot_regtargspairs++;

		$sr = $gs_size{$r}; $st = $gs_size{$t};
		$gs_term{$r} = $line[2]; $gs_term{$t} = $line[4];
		
		if($line[0]=~/R-C/) {
			$reg_gs{$r}++; $core_gs{$t}++;

			my $c = 0;		# No. common genes
			my $num_genepairs = 0;	# No. gene pairs across reg-tar terms
			foreach my $g1 (keys %{$gs_gene{$r}}) {
				foreach my $g2 (keys %{$gs_gene{$t}}) {
					if($g1 eq $g2) { $c++; }
					else {
						my $edge = join '__', sort($g1, $g2);
						$regtar_genepair{$edge}++;
						$num_genepairs++;
					}
				}
			}

			if($num_genepairs >= $minedges) {
				$num_regtargspairs_withenoughgenepairs++;
				$regtar_gspair{$r.'_'.$t}=$line[0];
				print ART "$r\t$sr\t$line[2]\t$t\t$st\t$line[4]\t$c\t$num_genepairs\n";
			}
		}
	}
}

my @reg_gsarray = (keys %reg_gs); my @core_gsarray = (keys %core_gs);
my ($num_reggs, $num_coregs) = (0, 0);
my %reg_genes = (); my %core_genes = ();
foreach my $r (keys %reg_gs) { $num_reggs++; foreach my $g (keys %{$gs_gene{$r}}) { $reg_genes{$g}++; } }
foreach my $c (keys %core_gs) { $num_coregs++; foreach my $g (keys %{$gs_gene{$c}}) { $core_genes{$g}++; } }
my $num_reggenes = scalar(keys %reg_genes); my $num_coregenes = scalar(keys %core_genes);

my %neg_regtar_genepair = ();
foreach my $r (keys %reg_gs) {
	foreach my $t (keys %core_gs) {
		unless(exists $regtar_gspair{$r.'_'.$t}) {
			foreach my $g1 (keys %{$gs_gene{$r}}) {
				foreach my $g2 (keys %{$gs_gene{$t}}) {
					my $edge = join '__', sort($g1, $g2);
					unless(($g1 eq $g2)||(exists $coann_genepair{$edge})||(exists $regtar_genepair{$edge})) {
						$neg_regtar_genepair{$edge}++;
					}
				}
			}
		}
	}
}

my @genes; my $num_posgenepairs = 0;
foreach my $edge (keys %regtar_genepair) {
	$num_posgenepairs++;
	@genes = split("__", $edge);
	print GS "$genes[0]\t$genes[1]\t1\n";
}

my $num_neggenepairs = 0;
foreach my $edge (keys %neg_regtar_genepair) {
	$num_neggenepairs++;
	@genes = split("__", $edge);
	print GS "$genes[0]\t$genes[1]\t0\n";
}

close GS; close ART;

print "\nTot. no. reg-tar gs pairs: $tot_regtargspairs\nNo. reg-tar pairs with enough edges: $num_regtargspairs_withenoughgenepairs\n";
print "\nNo. regulatory gs: $num_reggs\nNo. regulatory genes: $num_reggenes\n";
print "\nNo. core gs: $num_coregs\nNo. core genes: $num_coregenes\n";
print "\nNo. positive genepairs: $num_posgenepairs\nNo. negative genepairs: $num_neggenepairs\n\n";

# sub get_ovlp {}

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
