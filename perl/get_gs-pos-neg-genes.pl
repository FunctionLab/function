#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

sub parse_gmt;

my($igmt, $ianc, $igslist) = @ARGV;
my (@p, @ref);

@ref = parse_gmt($igmt);
my %gs_posg = %{$ref[0]};
my %gs_desc = %{$ref[1]};
my %gs_size = %{$ref[2]};

my %gs_bg = ();
open BGS, "$igslist" or die "Can't open $igslist!";
while(<BGS>) {
    if($_ =~ /^#/) { next; }
    chomp($_); $gs_bg{$_}++;
}
close BGS;

my %all_genes = ();
foreach (keys %gs_bg) {
    @all_genes{keys %{$gs_posg{$_}}} = values %{$gs_posg{$_}}; }

my %gs_anc = ();
open ANC, "$ianc" or die "Can't open $ianc!";
while(<ANC>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    unless((exists $gs_bg{$p[0]}) and (exists $gs_bg{$p[1]})) { next; }
    $gs_anc{$p[0]}{$p[1]} = $p[2];
}
close ANC;

my (%gs_negg, $par, $file);
foreach my $gs (keys %gs_bg) {
    %gs_negg = ();
    @gs_negg{keys %all_genes} = values %all_genes;

    foreach my $g (keys %all_genes) {
        $par = 0;

        if(exists $gs_posg{$gs}{$g}) { $par = 1; }
        else {
            foreach my $ags (keys %{$gs_anc{$gs}}) {
                if(exists $gs_posg{$ags}{$g}) {
                    $par = 1; last; } } }

        if($par == 0) { $gs_negg{$g}++; }
    }

    $file = $gs.'-gstd.txt'; print "$file\n";
    open GSTD, ">$file";
    foreach my $g (keys %{$gs_posg{$gs}}) { print GSTD "$g\t1\n"; }
    foreach my $g (keys %gs_negg) { print GSTD "$g\t-1\n"; }
    close GSTD;
}

# Assigns genes, description & size to gs
sub parse_gmt {
    my $gmt = shift;

    my %genes = (); my %desc = (); my %size = ();
    my (@q, $gs);

    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        @q = split '\t', $_;
        $gs = shift(@q);
        $q[0] =~ s/ \([0-9]*\)//g;
        $desc{$gs} = shift(@q);

        foreach my $g (@q) {
            $genes{$gs}{$g}++; }

        $size{$gs} = scalar keys %{$genes{$gs}};
    }
    close GMT;

    return (\%genes, \%desc, \%size);
}


