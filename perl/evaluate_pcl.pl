#!/usr/bin/perl

# ====================================================
# Name : evaluate_pcl.pl
# Purpose : Convert PCL to DAB and compare to gold-std
# Created : 26-09-2012
# Last Modified : Mon 08 Oct 2012 12:02:46 PM EDT
# Author(s) : Arjun Krishnan
# ====================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my ($help, $ipcl, $igs, $itag, $otab);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 4);
GetOptions( 'help' => \$help,
          'ipcl=s' => \$ipcl,
           'igs=s' => \$igs,
          'itag=s' => \$itag) or pod2usage( -exitstatus => 2, -verbose => 2 );
      # 'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

# my @measure = qw(pearson euclidean kendalls spearman pearnorm hypergeom innerprod bininnerprod mutinfo dcor sdcor);
my @measure = qw(pearson euclidean kendalls spearman pearnorm hypergeom innerprod bininnerprod mutinfo);
# my @measure = qw(pearson euclidean);
my $distancer = '/Genomics/fgrid/function/sleipnir/tools/Distancer/Distancer';
my $dchecker = '/Genomics/fgrid/function/sleipnir/tools/DChecker/DChecker';

if($igs =~ /\.dat$/) {
    my $dat = $igs; $igs =~ s/\.dat$/\.dab/g;
    `Dat2Dab -i $dat -Z -o $igs`; }

($otab = $ipcl) =~ s/\.pcl/\.$itag\.eval/g;
open TAB, ">$otab"; print TAB "Cut";

my ($time, @p, $dab, $eval);
my %prc = ();
foreach my $mes (@measure) {
    $time = runtime(); print "\n$mes";
    print TAB "\t$mes.RC\t$mes.PR";

    ($dab = $ipcl) =~ s/\.pcl/\.$mes\.dab/g;
    $time = runtime(); print "\n\t$time: $dab";
    unless(-e $dab) {
        `Distancer -i $ipcl -o $dab -d $mes`; }

    ($eval = $dab) =~ s/\.dab/\.$itag\.eval/g; $eval =~ s/^.*\///g;
    $time = runtime(); print "\n\t$time: $eval";
    unless(-e $eval) {
        `DChecker -w $igs -i $dab > $eval`; }

    open EVAL, "$eval" or die "Can't open $eval!";
    while (<EVAL>) {
        if($_ =~ /^#/) { next; }
        if($_ =~ /Cut/) { next; }
        chomp($_); @p = split '\t', $_;
        push(@{$prc{$p[0]}{$mes}}, $p[6], $p[7]);
        # print "\n$mes\t$eval\n\t@p\n\t${$prc{$p[0]}{$mes}}[0]\t${$prc{$p[0]}{$mes}}[1]\n";
        # exit;
    }
    close EVAL;
}
print TAB "\n";

my $ylim = 0;
foreach my $cut (sort {$b <=> $a} keys %prc) {
    print TAB "$cut";
    foreach my $mes (@measure) {
        print TAB "\t${$prc{$cut}{$mes}}[0]\t${$prc{$cut}{$mes}}[1]";
        if(${$prc{$cut}{$mes}}[1] > $ylim) {
            # print "\n$cut\t$mes\t${$prc{$cut}{$mes}}[1]\n"; exit;
            $ylim = ${$prc{$cut}{$mes}}[1]; }
    }
    print TAB "\n";
}

close TAB;

print "\t\n$ylim";
$ylim = sprintf("%.1f", ($ylim + 0.05)); print " --> $ylim\n";
`/home/arjunk/src/r/plot_prc_table.R $otab $otab.pdf $ylim`;

$time = runtime(); print "\n$time: DONE\n\n";

__END__

=head1

Evaluate pair-wise data from input PCL in comparison to gold-std.

=head1 USAGE

evaluate_pcl.pl [--ipcl INPUT_PCL] [--igs GOLD-STD_DAB] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in a PCL, calculates a pair-wise dataset using each of the
distance measures: "pearson", "euclidean", "kendalls", "spearman", "pearnorm",
"hypergeom", "innerprod", "bininnerprod", "mutinfo", "dcor", "sdcor". Then
conpares the resulting datasets to the gold-std and outputs precision-recall
points.

=head1 ARGUMENTS

=over 12

=item C<--ipcl>

Input file in PCL format.

=item C<--igs>

Input gold-std. If in DAB format, then its should contain all positive and
negative pairs. If in DAT format, its converted to DAB and negatives are
filled-in.

=item C<--itag>

Name tag denoting the gold-std that is used to define the output table
containing precision-recall points for each distance measure.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Sep 26

=cut

