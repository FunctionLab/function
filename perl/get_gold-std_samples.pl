#!/usr/bin/perl

# ================================================
# Name : get_gold-std_samples.pl
# Purpose : Get samples of the gold-std for evaluation
# Created : 01-02-2013
# Last Modified : Thu 14 Feb 2013 12:42:07 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my ($help, $istd); my $ifrac = 0.25; my $iprior = 0.05; my $inum = 10;
pod2usage( -exitstatus => 2, -verbose => 2) if (@ARGV < 1);
GetOptions( 'help' => \$help,
          'istd=s' => \$istd,
         'ifrac=f' => \$ifrac,
        'iprior=f' => \$iprior,
          'inum=i' => \$inum) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my ($time, @p);

unless($istd =~ /\.dab$/) { die "Gold-std should be in DAB format"; }

(my $pos_dat = $istd) =~ s/\.dab$/\.pos\.dat/g;
(my $neg_dab = $istd) =~ s/\.dab$/\.neg\.dab/g;

unless(-e $pos_dat) { `Dat2Dab -i $istd -c 1 -o $pos_dat`; }
unless(-e $neg_dab) { `Filterer -i $istd -o $neg_dab i=0`; }

chomp( my $tot_genes = `Dat2Dab -i $istd -E | wc -l` );
chomp( my $tot_pose = `wc -l $pos_dat` );
@p = split ' ', $tot_pose; $tot_pose = $p[0];

my $nnege = ((1-$iprior)/$iprior)*$ifrac*$tot_pose;

my ($sample_std, $pos_sdat, $fr, $neg_sdat);
for(my $i=1; $i<=$inum; $i++) {
    ($sample_std = $istd) =~ s/dab$//g; $sample_std =~ s/^.*\///g;
    $sample_std .= 'sample'.sprintf("%02s", $i).'.dat';

    $time = runtime(); print "\n$time: $sample_std";

    ($pos_sdat = $sample_std) =~ s/\.dat$/\.pos\.dat/g;
    ($neg_sdat = $sample_std) =~ s/\.dat$/\.neg\.dat/g;

    `Dat2Dab -i $pos_dat -u $ifrac -o $pos_sdat`;
    `Dat2Dab -i $pos_sdat -E > tmp.genes`;

    `Dat2Dab -i $neg_dab -g tmp.genes -o tmp.neg.dab`;
    chomp( $fr = `Dat2Dab -i tmp.neg.dab -P | head -1` );
    if($fr < $nnege) { die "$i not enough -ves: $fr < $nnege"; }
    $fr = sprintf("%.6f", ($nnege/$fr));
    `Dat2Dab -i $neg_dab -g tmp.genes -u $fr -o $neg_sdat`;

    `cat $pos_sdat $neg_sdat > $sample_std`;

    `rm -f $pos_sdat $neg_sdat tmp.genes tmp.neg.dab`;
}

$time = runtime(); print "\n$time: DONE\n\n";

__END__

=head1

Generate samples of gold-std for evaluation.

=head1 USAGE

get_gold-std_samples.pl [--istd GOLD-STD_DAB] [--ifrac FRAC_POS_DOUBLE] [--iprior
PRIOR_DOUBLE] [--inum NUM_SAMPLES_INT] [--help]

=head1 DESCRIPTION

This script takes in a binary gold-std and samples subsets for evaluation.

=head1 ARGUMENTS

=over 12

=item C<--istd>

Input gold-std file in DAB format.

=item C<--ifrac>

Fraction of positives to sample for a given subset.

=item C<--iprior>

Ratio of +ves to the total no. of instances (+ves plus -ves).

=item C<--inum>

No. of samples to generate.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 Feb 01

=cut

