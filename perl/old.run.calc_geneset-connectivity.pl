#!/usr/bin/perl

# ================================================
# Name : run_geneset-connectivity.pl
# Purpose : Run script to get the connectivity of network genes to multiple genesets
# Created : 01-05-2014
# Last Modified : Mon 23 Mar 2015 05:56:42 PM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

sub index_genes;

my $dat2dab = '/Genomics/fgrid/function/sleipnir-build/bin/Dat2Dab';
#my $conn = '/Genomics/ogtr04/arjunk/src/perl/calc_geneset-connectivity.pl';
#my $summ = '/Genomics/ogtr04/arjunk/src/perl/summarize-genes_multi-gs.pl';
my $conn = '/Genomics/ogtr04/arjunk/bin/calc_geneset-connectivity';
#my $summ = '/Genomics/ogtr04/arjunk/bin/summarize-genes_multi-gs';
my $summ = '/Genomics/ogtr04/arjunk/bin/calc_geneset-significance';

my ($help, @idaba, $igenes, $igmt1, $igmt2, $iing, $iexg, $iexe, $ibgcor, $iwseed, $iself, $inosumm, $ipara);
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
         'igenes=s' => \$igenes,
          'igmt1=s' => \$igmt1,
          'igmt2=s' => \$igmt2,
           'iing=s' => \$iing,
           'iexg=s' => \$iexg,
           'iexe=s' => \$iexe,
         'iming1=i' => \$iming1,
         'imaxg1=i' => \$imaxg1,
         'iming2=i' => \$iming2,
         'imaxg2=i' => \$imaxg2,
            'iself' => \$iself,
          'iuser=s' => \$iuser,
           'iwseed' => \$iwseed,
          'inosumm' => \$inosumm,
            'ipara' => \$ipara,
         'iqueue=s' => \$iqueue,
           'ibgcor' => \$ibgcor    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

unless($igmt2) {
    $igmt2 = $igmt1; }

if($iself) {
    $summ = '/Genomics/ogtr04/arjunk/src/perl/summarize-genes_single-gs.pl';
    $iqueue = '1hr'; }


my (%inclg, %exclg);
if($iing) {
    %inclg = index_genes($iing); }
if($iexg) {
    %exclg = index_genes($iexg); }


unless($igenes) {
    my $idab = $idaba[0];
    ($igenes = $idab) =~ s/\.da[tb]$/\.genes/g; $igenes =~ s/^.*\///g;
    unless(-e $igenes) {
        `$dat2dab -i $idab -E > $igenes`; } }
my %netg = index_genes($igenes);

(my $idegree = $idaba[0]) =~ s/\.da[tb]$/\.degree/g; $idegree =~ s/^.*\///g;
if($ibgcor) {
    if($#idaba > 0) {
        die "Currently can't do background correction for multiple networks. Run one at a time\n" }
    my $idab = $idaba[0];
    unless(-e $idegree) {
        `$dat2dab -i $idab -P | tail -n +3 > $idegree`; } }


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
        my $iconn = $tag.'.conn';
        my $osumm = $tag.'.summ';

        my @time = localtime(time);
        my $conn_job = 'conn.'.$tag.'.'. (join '',@time[0..3]);
        my $summ_job = 'summ.'.$tag.'.'. (join '',@time[0..3]);


        #perl ~/src/perl/calc_geneset-connectivity.pl --idab
        #../../../../net-cluster/v1/brain.degnorm.dab --igenew
        #../../brainspan/autism_prediction_rscores.exp100.txt --iseed
        #MP__0006018.genes --igmt
        #../../../../current/gene-sets/brainspan.st-zge2.gmt --iming 10 --imaxg
        #100000

        #perl ~/src/perl/calc_geneset-significance.pl --igene_scores
        #MP__0006018.brain.degnorm_genes.conn --igs_scores
        #MP__0006018.brain.degnorm_gs.conn --isize 202 --otab temp


        if($ibgcor) {
            unless(-e $iconn) {
                if($iexe) {
                    `qsub -m n -N $conn_job -l $iqueue -cwd "$conn --idab $idab --iseed $iseed --igdeg $idegree --iexe $iexe"`; }
                else {
                    `qsub -m n -N $conn_job -l $iqueue -cwd "$conn --idab $idab --iseed $iseed --igdeg $idegree"`; } }
            unless($inosumm) {
                unless(-e $osumm) {
                    `qsub -m n -hold_jid $conn_job -N $summ_job -l $iqueue -cwd "$summ --itab $iconn --igmt $igmt2 --iming $iming2 --imaxg $imaxg2 --otab $osumm"`; } } }
        else {
            unless(-e $iconn) {
                my $iextra_arg = '';

                if($iexe) {
                    $iextra_arg .= "--iexe $iexe"; }

                if($iwseed) {
                    $iextra_arg .= " --iwseed"; }

                `qsub -m n -N $conn_job -l $iqueue -cwd "$conn --idab $idab --iseed $iseed --igenes $igenes $iextra_arg"`; }

            unless($inosumm) {
                unless(-e $osumm) {
                    my $iextra_arg = '';

                    if($iself) {
                        if($iing) {
                            $iextra_arg .= "--iing $iing "; }
                        $iextra_arg .= "--igs $iseed --otxt $osumm"; }
                    else {
                        $iextra_arg .= "--igmt $igmt2 --iming $iming2 --imaxg $imaxg2 --otab $osumm"; }

                    if($ipara) {
                        $iextra_arg .= " --ipara"; }

                    #print "qsub -m n -hold_jid $conn_job -N $summ_job -l $iqueue -cwd \"$summ --itab $iconn $iextra_arg\""; exit;
                    `qsub -m n -hold_jid $conn_job -N $summ_job -l $iqueue -cwd "$summ --itab $iconn $iextra_arg"`;
                }
            }
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
[--ibgcor BACKGROUND_CORR] [--help]

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

=item C<--ibgcor>

(Optional) Option to calculate connectivity by correcting for the background of
the network gene. Default off.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

