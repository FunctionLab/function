#!/usr/bin/perl

# ================================================
# Name : run.map_list-to-list.pl
# Purpose : Run map_list-to-list.pl in parallel in cetus
# Created : 19-02-2015
# Last Modified : Sat 19 Mar 2016 11:15:48 AM EDT
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;

my $imap = '/Genomics/ogtr04/arjunk/bin/map_list-to-list';

my $injobs = 800;
my $iqueue = '1hr';
my $icut = 0.9;
my $icomb = 'max';


my $ilist1 = $ARGV[0];
my $ilist2 = $ARGV[1];


my %l1t = read_list($ilist1);


my $jcount; my $n = 0;
foreach my $l1 (keys %l1t) {
    (my $ilf1 = $l1); $ilf1 =~ s/:/__/g;

    chomp($jcount = `qstat -u arjunk | wc -l`);

    while($jcount > $injobs) {
        print "\t$jcount: sleeping for 2min ...\n";
        sleep 120;
        chomp($jcount = `qstat -u arjunk | wc -l`); }

    my $omap = $ilf1.'.map';
    my @time = localtime(time);
    my $jobid = $omap.'.'.(join '',@time[0..3]);

    $n++; print "$n\t$l1\n";
    unless(-e $omap) {
        #`qsub -m n -N $jobid -l $iqueue -cwd "$imap --ilist1 $ilf1 --ilist2 $ilist2 --icut $icut --icomb $icomb --otab $omap"`; }
        # `qsub -m n -N $jobid -l $iqueue -cwd "$imap --ilist1 $ilf1 --ilist2 $ilist2 --ionly1to2 --otab $omap"`; }
        `qsub -m n -N $jobid -l $iqueue -cwd "$imap --ilist1 $ilf1 --ilist2 $ilist2 --ionly1to2 --inosplit --otab $omap"`; }
    #exit;
}


sub read_list {
    my $ilist = shift;

    my %terms = ();
    open TR, "$ilist" or die "Can't open $ilist!";
    while (<TR>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;

        $terms{$p[0]}{$_}++; }
    close TR;

    foreach (keys %terms) {
        (my $ilf = $_); $ilf =~ s/:/__/g;
        unless(-e $ilf) {
            open LF, ">$ilf";
            foreach my $line (keys %{$terms{$_}}) {
                print LF "$line\n"; }
            close LF; } }

    return %terms; }

