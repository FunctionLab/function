#!/usr/bin/perl

# ================================================
# Name : calc_intra-geneset-connectivity_part2.pl
# Purpose : Compare real and random genesets to calculate significance
# Created : 04-01-2015
# Last Modified : Sun 11 Jan 2015 11:47:52 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($help, $idesc, $igdir, $ihdir, $otab);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 3);
GetOptions( 'help' => \$help,
         'idesc=s' => \$idesc,
         'igdir=s' => \$igdir,
         'ihdir=s' => \$ihdir,
          'otab=s' => \$otab    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %rgs_cliq = my %rgs_num = my %rgs_mean = my %rgs_sd = ();
chomp(my @randh = `ls $ihdir/s[0-9]*.r[0-9]*`);
foreach (@randh) {
    (my $rgs = $_) =~ s/^.*\///g;
    (my $n = $rgs) =~ s/\.r[0-9]*$//g; $n =~ s/^s//g;

    chomp(my $line = `grep s$n $_`);

    my @p = split '\t', $line;
    my $cliq = $p[$#p-2];
    push(@{$rgs_cliq{$n}}, $cliq);

    unless(exists $rgs_num{$n}) {
        $rgs_num{$n} = $rgs_mean{$n} = $rgs_sd{$n} = 0; }

    my $oldm = $rgs_mean{$n};

    $rgs_num{$n}++;
    $rgs_mean{$n} += ($cliq - $oldm)/$rgs_num{$n};
    $rgs_sd{$n} += ($cliq - $oldm)*($cliq - $rgs_mean{$n}); }


foreach my $n (sort keys %rgs_cliq) {
    my @a = sort {$b <=> $a} @{$rgs_cliq{$n}};
    $rgs_cliq{$n} = \@a;

    $rgs_sd{$n} = sqrt($rgs_sd{$n} / $rgs_num{$n}); }


my %gs_desc = my %gs_size = my %gs_cliq = my %gs_zscore = my %gs_pval = ();
open DESC, "$idesc" or die "Can't open $idesc!";
while (<DESC>) {
    if($_ =~ /^#/) { next; }
    chomp($_);
    
    my ($gs, $desc) = split '\t', $_;
    (my $gsid = $gs) =~ s/:/__/g;
    ($gs_desc{$gs} = $desc) =~ s/ \([0-9]*\)$//g;
    print "\n$gs\t$desc";

    my $gfile = $igdir.'/'.$gsid;
    unless(-e $gfile) { next; }
    open GEN, "$gfile";
    chomp(my @genes = <GEN>); close GEN;
    my $n = scalar @genes;
    $gs_size{$gs} = $n;

    my $hfile = $ihdir.'/'.$gsid;
    unless(-e $hfile) { next; }
    chomp(my $line = `grep $gsid $hfile`);
    my @p = split '\t', $line;
    my $cliq = $p[$#p-2];
    $gs_cliq{$gs} = $cliq;

    $gs_zscore{$gs} = ($cliq - $rgs_mean{$n}) / $rgs_sd{$n};
    $gs_pval{$gs} = emp_pval($cliq, $rgs_cliq{$n}); }
close DESC;
print "\n\n";

my %gs_esq = get_esfdr(\%gs_pval);


open TAB, ">$otab";
print TAB "#gs.id\tgs.size\tgs.desc\tgs.cliq\tgs.zscore\tgs.pval\tgs.esq\n";
foreach my $gs (sort keys %gs_pval) {
    print TAB "$gs\t$gs_size{$gs}\t$gs_desc{$gs}\t$gs_cliq{$gs}\t";
    printf TAB "%.6f\t%.6g\t%.6g\n", $gs_zscore{$gs}, $gs_pval{$gs}, $gs_esq{$gs}; }
close TAB;



# Get P-values based on empirical distribution
sub emp_pval {
    my $score = shift;
    my $sorts = shift;

    my $exc = 1; my $N = 1;
    foreach my $s (@$sorts) {
        $N++;
        if($s >= $score) {
            $exc++; } }

    return ($exc / $N); }


# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $pref = shift;
    
    my $minp = 1; my $ntests = 0;
    foreach (keys %{$pref}) {
        $ntests++;
        if(($pref->{$_} != 0) 
                and ($pref->{$_} < $minp)) {
            $minp = $pref->{$_}; } }

    my %esfdr = (); my $rank = $ntests;
    my $prev_qval = my $curr_qval = 1;

    foreach (sort {$pref->{$b} <=> $pref->{$a}} keys %{$pref}) {
        my $pvalue = $pref->{$_};

        if($pvalue == 0) {
            $curr_qval = ($minp * $ntests) / ($rank * 10); }
        else {
            $curr_qval = ($pvalue * $ntests) / $rank; }

        my $qvalue = 1;
        if($prev_qval < $curr_qval) {
            $qvalue = $prev_qval; }
        else {
            $qvalue = $curr_qval;
            $prev_qval = $curr_qval; }

        $esfdr{$_} = (-1)*log($qvalue)/log(10);
        if($esfdr{$_} < 0) { $esfdr{$_} = 0; }

        $rank--; }

    return %esfdr; }




__END__

=head1

Calculate significance of intra-geneset connectivity.

=head1 USAGE

calc_intra-geneset-connectivity_part2.pl [--idesc GENESET_DESC] [--igdir DIR_GENELISTS] [--ihdir DIR_HUBBER] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script calculates the significance of real genesets by comparing their
connectivity to those of random genesets of the same size.

=head1 ARGUMENTS

=over 12

=item C<--idesc>

File containing table of real genests along with their descriptions,
one-per-line: <gs.id> <desc>

=item C<--igdir>

Directory containing lists of genesets.

=item C<--ihdir>

Directory containing hubber result for each geneset.

=item C<--otab>

Output file containing geneset statistics.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2015 Jan 02

=cut

