#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my($in_list, $in_col, $in_idmap, $out_list) = @ARGV;
open FH, "$in_list" or die "Can't open file: $in_list!";
    chomp(my @f=<FH>); close FH;
open GH, "$in_idmap" or die "Can't open file: $in_idmap!";
    chomp(my @g=<GH>); close GH;
open HH, ">$out_list";

my (@p, $oid, $nid); my %map = ();
foreach (@g) {
    if($_ =~ /^#/) { next; }
    $_ =~ s///g;
    @p = split '\t', $_;
    $oid = shift @p;
    foreach my $nid (@p) { push(@{$map{$oid}}, $nid); }
}

my $tot_oid = 0; my $num_oid_wmap = 0;
my %nid_list = ();
foreach (@f) {
    if($_ =~ /^#/) { next; }
    $_ =~ s///g;
    $tot_oid++;
    @p = split '\t', $_;

    if(exists $map{$p[$in_col]}) {
        $num_oid_wmap++;
        foreach my $nid (@{$map{$p[$in_col]}}) {
            $nid_list{$nid} = $p[2];
        }
    }
}

my $num_nid = 0;
foreach my $nid (keys %nid_list) {
    $num_nid++;
    print HH "$nid\t$nid_list{$nid}\n";
}

print "\nTot. no. original ids: $tot_oid";
print "\nNo. orig ids with mapping: $num_oid_wmap";
print "\nNo. new ids: $num_nid\n\n";

close HH;

