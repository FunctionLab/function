#!/usr/bin/perl
use strict;
use warnings;

sub assign_col2to1 {
    my $array_ref = shift;
    my @p; my %map = (); my %count = ();

    foreach (@$array_ref) {
        if($_ =~ /^#/) { next; }
        @p = split '\t', $_;
        $map{$p[0]} = $p[1];
        $count{$p[0]}++;
    }

    foreach (keys %map) {
        if($count{$_} > 1) {
            delete $map{$_};
        }
    }

    return \%map;
}

my ($in_gmt, $out_gmt) = @ARGV;
my $in_idmap1 = './mappings/symbol_entrez_mapping_v2.txt';
my $in_idmap2 = './mappings/symbol-alias_entrez_mapping_v2.txt';

open FH, "$in_gmt" or die "Can't open $in_gmt";
    chomp(my @gmt=<FH>); close FH;
open P1H, "$in_idmap1" or die "Can't open $in_idmap1!";
    chomp(my @idmap1=<P1H>); close P1H;
open P2H, "$in_idmap2" or die "Can't open $in_idmap2!";
    chomp(my @idmap2=<P2H>); close P2H;
open GMT, ">$out_gmt";

my $hashref;
$hashref = assign_col2to1(\@idmap1);  my %id_map1 = %{$hashref};
$hashref = assign_col2to1(\@idmap2);  my %id_map2 = %{$hashref};

my (@p, $gs, $desc, $ngenes, %genes);
foreach (@gmt) {
    @p = split '\t', $_;
    $gs = shift(@p);
    $desc = shift(@p); $desc =~ s/ \([0-9]*\)$//g;
    
    %genes = ();
    foreach my $g (@p) {
        if(exists $id_map1{$g}) {
            $genes{$id_map1{$g}}++;
        }
        elsif(exists $id_map2{$g}) {
            $genes{$id_map2{$g}}++;
        }
    }

    $ngenes = scalar(keys %genes);
    print GMT "$gs\t$desc ($ngenes)";
    foreach my $g (keys %genes) {
        print GMT "\t$g";
    }
    print GMT "\n";
}

close GMT;

