#!/usr/bin/perl

# ================================================
# Name : run_geneset-connectivity.pl
# Purpose : Run script to get the connectivity of network genes to multiple genesets
# Created : 01-05-2014
# Last Modified : Thu 02 Apr 2015 05:47:32 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

sub index_genes;

my $dat2dab = '/Genomics/fgrid/function/sleipnir-build/bin/Dat2Dab';
my $conn = '/Genomics/ogtr04/arjunk/bin/calc_geneset-connectivity';
my $sigf = '/Genomics/ogtr04/arjunk/bin/calc_geneset-significance';
#my $conn = '/Genomics/ogtr04/arjunk/src/perl/calc_geneset-connectivity.pl';
#my $sigf = '/Genomics/ogtr04/arjunk/src/perl/calc_geneset-significance.pl';

my ($help, @idaba, $igenew, $igmt1, $igmt2,
    $iing, $iexg, $iexe, $ibgg, $ibgs, $iself,
    $inosigf, $ipara);
my $iming1 = my $iming2 = 5;
my $imaxg1 = my $imaxg2 = 100;

my $iuser = 'arjunk';

my $iqueue = '1hr';
#my $iqueue = '1day';

#my $injobs = 500;
my $injobs = 800;

pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(  'help' => \$help,
          'idaba=s' => \@idaba,
         'igenew=s' => \$igenew,
          'igmt1=s' => \$igmt1,
          'igmt2=s' => \$igmt2,
           'ibgg=s' => \$ibgg,
           'ibgs=s' => \$ibgs,
           'iing=s' => \$iing,
           'iexg=s' => \$iexg,
           'iexe=s' => \$iexe,
         'iming1=i' => \$iming1,
         'imaxg1=i' => \$imaxg1,
         'iming2=i' => \$iming2,
         'imaxg2=i' => \$imaxg2,
            'iself' => \$iself,
          'iuser=s' => \$iuser,
          'inosigf' => \$inosigf,
            'ipara' => \$ipara,
         'iqueue=s' => \$iqueue   ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

unless($igmt2) {
    $igmt2 = $igmt1; }

if($iself) {
    $sigf = '/Genomics/ogtr04/arjunk/src/perl/summarize-genes_single-gs.pl';
    $iqueue = '1hr'; }


my (%inclg, %exclg);
if($iing) {
    %inclg = index_genes($iing); }
if($iexg) {
    %exclg = index_genes($iexg); }


unless($igenew) {
    my $idab = $idaba[0];
    ($igenew = $idab) =~ s/\.da[tb]$/\.genes/g; $igenew =~ s/^.*\///g;
    unless(-e $igenew) {
        `$dat2dab -i $idab -E > $igenew`; } }
my %netg = index_genes($igenew);


my %size_gs2 = ();
open GMT, "$igmt2" or die "Can't open $igmt2!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    (my $gs = $p[0]) =~ s/:/__/g;
    my $s = (scalar @p) - 2;
    if(($s < $iming2) or ($s > $imaxg2)) { next; }
    $size_gs2{$s}{$gs}++; }
close GMT;


my %size_bgs = ();
if($ibgs) {
    open BGS, "$ibgs" or die "Can't open $ibgs!";
    while (<BGS>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $size_bgs{$p[0]} = $p[1]; }
    close BGS; }


print "\n";
my $ngs = 0; my $jcount;
open GMT, "$igmt1" or die "Can't open $igmt1!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    (my $gs = $p[0]) =~ s/:/__/g;
    shift @p;

    my %gsg = ();
    foreach my $g (@p) {
        unless(exists $netg{$g}) { next; }
        if($iing) { unless(exists $inclg{$g}) { next; } }
        if($iexg) { if(exists $exclg{$g}) { next; } }
        $gsg{$g}++; }

    my $size = scalar keys %gsg;
    if(($size < $iming1) or ($size > $imaxg1)) { next; }

    $ngs++;
    print "$ngs\t$gs\t$size\n";

    my $iseed = $gs.'.genes';
    unless(-e $iseed) {
        open SET, ">$iseed";
        foreach my $g (sort keys %gsg) {
            print SET "$g\n"; }
        close SET; }

    foreach my $idab (@idaba) {
        chomp($jcount = `qstat -u $iuser | wc -l`);
        while($jcount > $injobs) {
            print "\t$jcount: sleeping for 1min ...\n";
            sleep 60;
            chomp($jcount = `qstat -u $iuser | wc -l`); }

        (my $tag = $idab) =~ s/^.*\///g; $tag =~ s/\.da[tb]$//g;
        $tag = $gs.'.'.$tag;
        my $iconn_genes = $tag.'_genes.conn';
        my $iconn_gs = $tag.'_gs.conn';

        my @time = localtime(time);
        my $conn_job = 'conn.'.$tag.'.'. (join '',@time[0..3]);

        unless(-e $iconn_genes) {
            my $iextra_arg = '';

            if($iexe) {
                $iextra_arg .= "--iexe $iexe"; }

            if($ibgg) {
                $iextra_arg .= " --ibgg $ibgg"; }

            #print "\nqsub -m n -N $conn_job -l $iqueue -cwd \"$conn --idab $idab --iseed $iseed --igenew $igenew --igmt $igmt2 --iming $iming2 --imaxg $imaxg2 $iextra_arg\""; }
            `qsub -m n -N $conn_job -l $iqueue -cwd "$conn --idab $idab --iseed $iseed --igenew $igenew --igmt $igmt2 --iming $iming2 --imaxg $imaxg2 $iextra_arg"`; }

        if($inosigf) { next; }

        foreach my $s (keys %size_gs2) {
            my $osigf = $tag.'.'.$s.'.sigf';
            if(-e $osigf) { next; }

            chomp($jcount = `qstat -u $iuser | wc -l`);
            while($jcount > $injobs) {
                print "\t$jcount: sleeping for 2min ...\n";
                sleep 120;
                chomp($jcount = `qstat -u $iuser | wc -l`); }

            my $sigf_job = join '.', ('sigf',$tag,$s,(join '',@time[0..3]));

            my $iextra_arg = '';

            if($iself) {
                if($iing) {
                    $iextra_arg .= "--iing $iing "; }
                $iextra_arg .= "--igs $iseed --otxt $osigf"; }
            else {
                $iextra_arg .= "--igene_scores $iconn_genes --igs_scores $iconn_gs --isize $s --otab $osigf"; }

            if($ibgs) {
                $iextra_arg .= " --ibgs $size_bgs{$s}"; }

            if($ipara) {
                $iextra_arg .= " --ipara"; }

            #print "\nqsub -m n -hold_jid $conn_job -N $sigf_job -l $iqueue -cwd \"$sigf $iextra_arg\""; exit;
            `qsub -m n -hold_jid $conn_job -N $sigf_job -l $iqueue -cwd "$sigf $iextra_arg"`;
        }
    }
}
close GMT;

print "\n$ngs genesets\n\n";



sub index_genes {
    my $ifile = shift;
    my %genes = ();

    open GEN, "$ifile" or die "Can't open $ifile!";
    while (<GEN>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $genes{$p[0]}++; }
    close GEN;

    return %genes; }


__END__

=head1

Run script to get connectivity of network genes to geneset in a collection.

=head1 USAGE

run_geneset-connectivity.pl [--idab NETWORK_DAB/DAT] [--igmt
GENESETS_GMT] [--iming MIN_GENES] [--imaxg MAX_GENES]
[--help]

=head1 DESCRIPTION

This script takes a network and a geneset collection as input and run another
script 'calc_geneset-connectivity.pl' to get, for each geneset,  a score for
each gene in the network that represents a connectivity score of that gene
to that geneset. Best followed by running the script
'collate_geneset-connectivity.pl' to put all the scores into a single matrix.

=head1 ARGUMENTS

=over 12

=item C<--idab>

Input network in DAB/DAT format.

=item C<--igmt1>

Input geneset collection in GMT format. This geneset is used to calcualte
connectivity of genes.

=item C<--igmt2>

(Optional) Second input geneset collection in GMT format. This geneset, if
provided, is used at the summarization stage. If not provided, genesets from the
first collection (--igmt1) will be used.

=item C<--iming>

(Optional) Min. no. genes in a geneset for analysis. Default 5.

=item C<--imaxg>

(Optional) Max. no. genes in a geneset for analysis. Default 200.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

