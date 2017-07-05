#!/usr/bin/perl

# ================================================
# Name : run.gsea_hypergeometric.pl
# Purpose : Run gsea_hypergeometric.pl in parallel
# Created : 20-02-2015
# Last Modified : Thu 09 Apr 2015 09:50:32 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

my $gsea = '/Genomics/ogtr04/arjunk/bin/gsea_hypergeometric';


my $igmt1 = $ARGV[0];
my $ibg1 = $ARGV[1];
my $igmt2 = $ARGV[2];

#(my $ibg1 = $igmt1) =~ s/\.gmt$/\.genes/g;


my $iming = 5;
my $imaxg1 = 500;
my $imaxg2 = 500;
my $imincom = 5;
#my $iqval = 0.01;
#my $ipval = 0.0001;
my $iovlp = 0.5;

my $injobs = 500;


my $jcount; my $n = 0;
open GMT, "$igmt1" or die "Can't open $igmt1!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    chomp($jcount = `qstat -u arjunk | wc -l`);
    while($jcount > $injobs) {
        print "\t$jcount: sleeping for 30sec ...\n";
        sleep 30;
        chomp($jcount = `qstat -u arjunk | wc -l`); }

    my $it1f = $p[0].'.gmt'; $it1f =~ s/:/__/g;
    my $ot1f = $p[0].'.ovp'; $ot1f =~ s/:/__/g;

    unless(-e $it1f) {
        open T1F, ">$it1f"; print T1F "$_\n"; close T1F; }

    my @time = localtime(time);
    my $jobid = $ot1f.'.'.(join '',@time[0..3]);
    $n++; print "$n\t$p[0]\t$p[1]\n";

    unless(-e $ot1f) {
        #`qsub -m n -N $jobid -l 1hr -cwd "$gsea --igmt1 $it1f --igmt2 $igmt2 --ibg1 $ibg1 --iming $iming --imaxg1 $imaxg1 --iqval $iqval --orich $ot1f"`;
        #`qsub -m n -N $jobid -l 1hr -cwd "$gsea --igmt1 $it1f --igmt2 $igmt2 --ibg1 $ibg1 --iming $iming --imaxg1 $imaxg1 --imaxg2 $imaxg2 --imincom $imincom --ipval $ipval --orich $ot1f"`;
        `qsub -m n -N $jobid -l 1hr -cwd "$gsea --igmt1 $it1f --igmt2 $igmt2 --ibg1 $ibg1 --iming $iming --imaxg1 $imaxg1 --imaxg2 $imaxg2 --iovlp $iovlp --orich $ot1f"`;
    }
}
close GMT;

