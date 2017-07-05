#!/usr/bin/perl

# ================================================
# Name : get_reciprocal-best-hits.pl
# Purpose : Process BLAST results to get reciprocal-best-hit pairs
# Created : 18-01-2013
# Last Modified : Fri 18 Jan 2013 04:39:34 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my($help, $im1, $im2, $otab); my $imes = 'iden';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
           'im1=s' => \$im1,
           'im2=s' => \$im2,
          'imes=s' => \$imes,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $iden, $bits); my %seq_len = ();

my %m1_best = (); my %m1_iden = (); my %m1_bits = ();

open IM1, "$im1" or die "Can't open $im1!";
while (<IM1>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $p[0] =~ s/\|.*$//g; $seq_len{$p[0]} = $p[1];
    $p[2] =~ s/\|.*$//g; $seq_len{$p[2]} = $p[3];

    $iden = (2*$p[4])/($p[1]+$p[3]);
    $bits = $p[$#p];

    if(exists $m1_best{$p[0]}) {
        if($imes eq 'iden') { if($iden <= $m1_iden{$p[0]}) { next; } }
        else { if($p[$#p] <= $m1_bits{$p[0]}) { next; } }
        $m1_best{$p[0]} = $p[2];
        $m1_iden{$p[0]} = $iden;
        $m1_bits{$p[0]} = $bits;
    }
    else {
        $m1_best{$p[0]} = $p[2];
        $m1_iden{$p[0]} = $iden;
        $m1_bits{$p[0]} = $bits; } }
close IM1;

my %m2_best = (); my %m2_iden = (); my %m2_bits = ();

open IM2, "$im2" or die "Can't open $im2!";
while (<IM2>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

    $p[0] =~ s/\|.*$//g; $seq_len{$p[0]} = $p[1];
    $p[2] =~ s/\|.*$//g; $seq_len{$p[2]} = $p[3];

    $iden = (2*$p[4])/($p[1]+$p[3]);
    $bits = $p[$#p];

    if(exists $m2_best{$p[0]}) {
        if($imes eq 'iden') { if($iden <= $m2_iden{$p[0]}) { next; } }
        else { if($p[$#p] <= $m2_bits{$p[0]}) { next; } }
        $m2_best{$p[0]} = $p[2];
        $m2_iden{$p[0]} = $iden;
        $m2_bits{$p[0]} = $bits;
    }
    else {
        $m2_best{$p[0]} = $p[2];
        $m2_iden{$p[0]} = $iden;
        $m2_bits{$p[0]} = $bits; } }
close IM2;

open TAB, ">$otab";
print TAB "#seq1.id\tseq1.len\tseq2.id\tseq2.len\tides\tbits\n";
# print TAB "#seq1.id\tseq1.len\tseq2.id\tseq2.len\tm1.ides\tm1.bits\tm2.ides\tm2.bits\n";
foreach my $s1 (keys %m1_best) {
    if($m2_best{$m1_best{$s1}} eq $s1) {
        my $s2 = $m1_best{$s1};

        print TAB "$s1\t$seq_len{$s1}\t$s2\t$seq_len{$s2}\t";
        print TAB sprintf("%.3f", $m1_iden{$s1}), "\t$m1_bits{$s1}\n"; } }
close TAB;


__END__

=head1

Calculates reciprocal best hits from two-way BLAST results.

=head1 USAGE

./get_reciprocal-best-hits.pl [--im1 IN_MATCH1] [--im2 IN_MATCH2] [--otab OUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in two BLAST results –– comparison between two sets of
sequences, bothways [one-vs-two and two-vs-one] –– and reports pairs of
sequences from the two sets of sequences that have best matches with each other.

	# Example BLAST run
	# =================
	
	If the two sequence datasets you wish to compare are seq1.fa and seq2.fa
	(sequence files in FASTA format), then, take the following steps:
	
	# Create BLAST databases
	makeblastdb -in seq1.fa -dbtype <nucl/prot> [do 'makeblastdb -help' for more info.]
	makeblastdb -in seq2.fa -dbtype <nucl/prot>
	
	# Run BLAST
	blastn -query seq1.fa -db seq2.fa -out 1-vs-2.txt -outfmt '6 qseqid qlen sseqid slen nident pident length mismatch positive gapopen qstart qend sstart send evalue bitscore'
	blastn -query seq2.fa -db seq1.fa -out 2-vs-1.txt -outfmt '6 qseqid qlen sseqid slen nident pident length mismatch positive gapopen qstart qend sstart send evalue bitscore'
	
	These two output files, 1-vs-2.txt and 2-vs-1.txt can then be processed with
	this script to generate reciprocally best matched sequences.

=head1 ARGUMENTS

=over 12

=item C<--im1>

First BLAST result in tabluar format. Generated with -outfmt '6 qseqid qlen
sseqid slen nident pident length mismatch positive gapopen qstart qend sstart
send evalue bitscore' option.

=item C<--im2>

Second BLAST result in tabluar format. Same format as --im1.

=item C<--imes>

(Optional) The measure – 'iden' or 'bits' – used to decide best match. If
'iden', then the number of identical matches relative to the lengths of the two
sequences is used as the measure. If 'bits', then the bti score from the BLAST
output is taken as the measure. Default 'iden'.

=item C<--otab>

Output table of reciprocal best matches.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 Jan 18

=cut

