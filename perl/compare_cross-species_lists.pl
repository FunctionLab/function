#!/usr/bin/perl
use strict;
use warnings;
use PDF;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

# Assigns family to genes
sub genes2fam {
    my $file = shift;
    my $s1 = shift;
    my $s2 = shift;

    open FH, "$file" or die "Can't open $file!";
        chomp(my @ortho=<FH>); close FH;

    my (@fl, $fam, $o, $g); my %gene_fam = ();
    foreach (@ortho) {
        unless(($_ =~ /$s1/) or ($_ =~ /$s2/)) { next; }
        @fl = split '\t', $_;
        $fam = shift(@fl);
        foreach my $h (@fl) {
            ($o, $g) = split /\|/, $h;
            $gene_fam{$g} = $fam;
        }
    }

    return %gene_fam;
}

# Assigns gene to family
sub list2fam {
    my $file = shift;
    my $famref = shift;
    open LH, "$file" or die "Ca't open $file!";
        chomp(my @list=<LH>); close LH;

    my %famlist = ();
    foreach (@list) {
        if(exists $famref->{$_}) {
            $famlist{$famref->{$_}}++;
        }
    }

    return %famlist;
}

my($in_list1, $in_list2, $in_fam, $in_sp, $out_table) = @ARGV;

my ($s1, $s2) = split ":", $in_sp;
my %gene_fam = genes2fam($in_fam, $s1, $s2);

my $univ = scalar values %gene_fam;

my %s1_famlist = list2fam($in_list1, \%gene_fam);
my %s2_famlist = list2fam($in_list2, \%gene_fam);

my ($common, $min, $union, $tag, $ovlp, $jac, $lor, $pvalue);

my %list_fam = ();
@list_fam{keys %s1_famlist} = values %s1_famlist;
@list_fam{keys %s2_famlist} = values %s2_famlist;
$union = scalar keys %list_fam;

my $s1_size = scalar keys %s1_famlist;
my $s2_size = scalar keys %s2_famlist;
$min = $s1_size; if($s2_size < $min) { $min = $s2_size; }

$common = 0;
foreach (keys %s1_famlist) {
    if(exists $s2_famlist{$_}) { $common++; }
}

$ovlp = ($common/$min);
$jac = ($common/$union);
$lor = ($common/$s1_size)/($s2_size/$univ);
unless($common == 0) { $lor = log($lor)/log(2); }
$pvalue = hypergeometric_tail($univ, $s1_size, $s2_size, $common);

print "\nCompared ...\n$in_list1 ($s1_size) ... and\n$in_list2 ($s2_size)\n";
print "\nUnion: $union\nCommon: $common\n";
print "\nOvlp: ", sprintf("%.3g", $ovlp), "\nJaccard: ", sprintf("%.3g", $jac), "\n";
print "LOR: ", sprintf("%.3g", $lor), "\nP-value: ", sprintf("%.3g", $pvalue), "\n\n";

