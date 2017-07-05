#!/usr/bin/perl

# ================================================
# Name : get_gold-std_for_list-of-terms.pl
# Purpose : Create a binary gold-std for list of terms
# Created : 13-01-2015
# Last Modified : Tue 13 Jan 2015 02:42:54 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

sub index_list;


my ($help, $isubgs, $igmt, $igenes, $odat);
my $iposgs = '/Genomics/ogtr04/arjunk/data/functional-annotations/go/goslim_2012.pos.txt';
my $ineggs = '/Genomics/ogtr04/arjunk/data/functional-annotations/go/goslim_2012.neg.txt';

pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 3);
GetOptions(     'help' => \$help,
            'isubgs=s' => \$isubgs,
            'ineggs=s' => \$ineggs,
              'igmt=s' => \$igmt,
            'igenes=s' => \$igenes,
              'odat=s' => \$odat      ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


print "\nindexing lists ...";
my %filt_genes = ();
if($igenes) {
    %filt_genes = index_list($igenes); }


my %sub_gs = index_list($isubgs);
my %pos_gs = index_list($iposgs);
my %neg_gs = index_list($ineggs);


print "\nindexing pos & slim edges ...";
my %std_edges = my %pos_genes = my %sub_genes = my %slim_gs_genes = ();
open GMT, "$igmt" or die "Can't open $igmt!";
while (<GMT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    my $gs = shift @p;
    my $desc = shift @p;

    if(exists $pos_gs{$gs}) {
        foreach my $g (@p) {
            if($igenes) {
                if(exists $filt_genes{$g}) {
                    $pos_genes{$g}++; } }
            else {
                $pos_genes{$g}++; } } }

    if(exists $sub_gs{$gs}) {
        print "\n\t+ $gs $desc";

        my @fgen = ();
        foreach my $g (@p) {
            if($igenes) {
                if(exists $filt_genes{$g}) {
                    $sub_genes{$g}++;
                    push(@fgen, $g); } }
            else {
                $sub_genes{$g}++;
                push(@fgen, $g); } }
        print " *", scalar @fgen; print "*";
       
        for(my $i=0; $i<$#fgen; $i++) {
            for(my $j=($i+1); $j<=$#fgen; $j++) {
                my $e = join '__', sort($fgen[$i], $fgen[$j]);
                $std_edges{$e} = 1; } } }

    elsif(exists $neg_gs{$gs}) {
        print "\n\t- $gs $desc $#p";

        foreach my $g (@p) {
            #for(my $i=0; $i<$#p; $i++) {
            #    for(my $j=($i+1); $j<=$#p; $j++) {
            #        my $e = join '__', sort($p[$i], $p[$j]);
            #        $slim_edges{$e}++; } }
            $slim_gs_genes{$gs}{$g}++; } } }
close GMT;


my %slim_ngenes = ();
foreach my $s (keys %slim_gs_genes) {
    $slim_ngenes{$s} = scalar keys %{$slim_gs_genes{$s}}; }

my @slima = sort {$slim_ngenes{$b} <=> $slim_ngenes{$a}} keys %slim_ngenes;


print "\nindexing neg edges ...";
my @subg = keys %sub_genes;
#for(my $i=0; $i<$#subg; $i++) {
#    for(my $j=($i+1); $j<=$#subg; $j++) {
        #if(exists $slim_edges{$e}) { next; }
foreach my $g (keys %pos_genes) {
    foreach my $h (@subg) {
        if($g eq $h) { next; }

        my $par = 0;
        foreach my $s (@slima) {
            if((exists $slim_gs_genes{$s}{$g})
                    and (exists $slim_gs_genes{$s}{$h})) {
                $par = 1; last; } }

        if($par == 1) { next; }

        my $e = join '__', sort($g, $h);
        $std_edges{$e} = 0; } }


my $npose = my $nnege = 0;

print "\nprinting gold-std ...";
open DAT, ">$odat";
foreach my $e (keys %std_edges) {
    my @p = split '__', $e;
    my $v = $std_edges{$e};

    print DAT "$p[0]\t$p[1]\t$v\n";

    if($v == 1) { $npose++; }
    elsif($v == 0) { $nnege++; } }
close DAT;

print "\n\nngene: ", scalar @subg;
print "\nnpose: $npose\nnnege: $nnege\n\n";



sub index_list {
    my $ilist = shift;
    my %index = ();

    open LIST, "$ilist" or die "Can't open $ilist!";
    while (<LIST>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $index{$p[0]}++; }
    close LIST;

    return %index; }


__END__

=head1

Get binary gold-std for list of terms.

=head1 USAGE

get_gold-std_for_list-of-terms.pl [--isubgs LIST_OF_POS-TERMS] [--ineggs
LIST_OF_NEG-TERMS] [--igmt GENESET_COLLECTION] [--igenes GENES_TO_INCLUDE]
[--odat OUTPUT_GOLD-STD] [--help]

=head1 DESCRIPTION

This script takes in a list of terms and gene-term associations, and outputs a
gold-std set of gene pairs, where gene pairs co-annotated to any of the input
terms are set to 1, and relevant gene pairs not co-annotated to a negative slim
are set to 0.

=head1 ARGUMENTS

=over 12

=item C<--isubgs>

Input list of terms to determine positive gene pairs in the gold-std: <id> <term>

=item C<--ineggs>

(Optional) Input list of terms to determine negative gene pairs in the gold-std:
<id> <term>. Default is the list of neagtive GO slim terms.

=item C<--igmt>

Input geneset collection in GMT format.

=item C<--igenes>

(Optional) List of genes to include in the gold-std.

=item C<--odat>

Output file containing the binary gold-std.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2015 Jan 13

=cut

