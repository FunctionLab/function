#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);

my ($help, $in_gstd, $in_glist_dir, $in_dat_dir, $out_ssv, $time);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
         'igstd=s' => \$in_gstd,
    'iglist-dir=s' => \$in_glist_dir,
      'idat-dir=s' => \$in_dat_dir,
          'ossv=s' => \$out_ssv) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

$in_glist_dir =~ s/\/$//g;
unless($in_dat_dir =~ /\/$/) { $in_dat_dir .= '/' }

open STD, "$in_gstd" or die "Can't open $in_gstd";
    chomp(my @std=<STD>); close STD;
chomp(my @dat_files = `ls $in_glist_dir/*.genes`);
open DAT, ">$out_ssv";

$time = runtime(); print "\n$time: Hashing dataset genes and gold-std edges...";
my @p; my %gstd_edges = ();
foreach my $e (@std) {
    @p = split '\t', $e;
    $e = join '__', sort($p[0], $p[1]);
    $gstd_edges{$e} = $p[2];
}

my $nume = scalar keys %gstd_edges; print "\n\tTot. edges: $nume";
my $nmin = ($nume/5); print "\n\tMin. edges to cover: $nmin\n";

my %dat_genes = ();
foreach my $dat (@dat_files) {
    open DG, "$dat" or die "Can't open $dat!";
        chomp(my @dg=<DG>); close DG;

    $dat =~ s/\.genes//g; $dat =~ s/^.*\///g;
    foreach my $g (@dg) {
        $dat_genes{$dat}{$g}++;
    }
}

my %sel_dat = (); my $edge_count; my @gstd_dat = (); my %sel_edges = ();
foreach my $dat (keys %dat_genes) {
    $edge_count = 0;
    foreach my $e (keys %gstd_edges) {
        @p = split '__', $e;
        if((exists $dat_genes{$dat}{$p[0]}) and (exists $dat_genes{$dat}{$p[1]})) {
            $edge_count++; $sel_edges{$e}++;
        }
    }

    if($edge_count < $nmin) { next; }
    push(@gstd_dat, $dat);
    $sel_dat{$dat}++;
}

print "\n\tNo. total datasets: ", scalar keys %dat_genes;
print "\n\tNo. selelcted datasets: ", scalar keys %sel_dat;

print "\n\n\tNo. total edges: ", scalar keys %gstd_edges;
print "\n\tNo. selected edges: ", scalar keys %sel_edges;

$time = runtime(); print "\n$time: Extracting data points ...";
my ($dab, $quant, $dab_w0s, $gsdat, $e);
my %gsdat_edges = (); my %max = ();
foreach my $dat (keys %sel_dat) {
    $dab = $in_dat_dir.$dat.'.qdab';
    # ($quant = $dab) =~ s/\.qdab/\.quant/g;
    # $dab_w0s = $dat.'.w0s.qdab';
    $gsdat = $dat.'-gs.dat';

    $time = runtime(); print "\n\t$time: $dat";
    # `Dat2Dab -i $dab -q $quant -Z -o $dab_w0s`;
    # `Dat2Dab -i $dab_w0s -e $in_gstd -o $gsdat`;
    `Dat2Dab -i $dab -e $in_gstd -o $gsdat`;
    
    $max{$dat} = 0;
    open DH, "$gsdat"; chomp(my @gsdat=<DH>); close DH;
    foreach my $pair (@gsdat) {
        @p = split '\t', $pair;
        $e = join '__', sort($p[0], $p[1]);
        $gsdat_edges{$dat}{$e} = $p[2];
        if(abs($p[2]) > $max{$dat}) { $max{$dat} = abs($p[2]); }
    }
    `rm -f $dab_w0s`;
    `rm -f $gsdat`;
}

$time = runtime(); print "\n$time: Printing ssv's ...\n";
my ($dat, @dat_array); $nume = 0;

foreach my $e (keys %gstd_edges) {
    unless(exists $sel_edges{$e}) { next; }
    $nume++; if(($nume % 10000) == 0) { print "\t$time: $nume edges\n"; }
    print DAT "$gstd_edges{$e}";

    for(my $j=0; $j<=$#gstd_dat; $j++) {
        $dat = $gstd_dat[$j];
        unless(exists $gsdat_edges{$dat}{$e}) { next; }
        print DAT " ", ($j+1), ":", sprintf("%.3f",
            ($gsdat_edges{$dat}{$e}/$max{$dat}));
    }
    print DAT "\n";
}

close DAT;

$time = runtime();
print "$time: DONE\n\n";


__END__

=head1 USAGE

./get_data_for_gold-std.pl [--igstd GOLD-STD_DAT] [--iglist-dir GENE-LIST_DIR]
[--idat-dir DATASETS_DIR] [--ossv OUTPUT_SSV] [--help]

=head1 DESCRIPTION

This script takes in a gold-std DAT file and appends edge-level data from a
collection of DAB files. The output format is an SSV, which can be used directly
for classification or regression.

=head1 ARGUMENTS

=over 12

=item C<--igstd>

Gold-std edges. DAT format.

=item C<--iglist-dir>

Directory containing the list of genes in each DAB in the dataset collection.

=item C<--idat-dir>

Directory containing the dataset collection (DABs).

=item C<--ossv>

Output file. SSV format.

=item C<--help>

prints this documentation

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Apr 24

=cut

