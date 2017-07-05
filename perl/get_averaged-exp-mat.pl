#!/usr/bin/perl

# ================================================
# Name : get_averaged-exp-mat.pl
# Purpose : Get expression matrix with values averaged across replicates
# Created : 27-09-2013
# Last Modified : Mon 08 Feb 2016 05:10:23 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

sub get_summary;


my ($help, $imat, $isam, $isiz, $ifil, $igen, $icut, $ogrp, $omat, $ovar); my $iskip = 0; my $istat = 'mean';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'isam=s' => \$isam,
          'isiz=s' => \$isiz,
          'ifil=s' => \$ifil,
          'igen=s' => \$igen,
          'icut=s' => \$icut,
         'iskip=i' => \$iskip,
         'istat=s' => \$istat,
          'ogrp=s' => \$ogrp,
          'ovar=s' => \$ovar,
          'omat=s' => \$omat ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %incl_samp = ();
if($ifil) {
    open FIL, "$ifil" or die "Can't open $ifil!";
    while (<FIL>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $incl_samp{$_}++; }
    close FIL; }


my %sample_siz = ();
if($isiz) {
    open SIZ, "$isiz" or die "Can't open $isiz!";
    while (<SIZ>) {
        if($_ =~ /^#/) { next; }
        if($_ =~ /^group/) { next; }
        chomp($_); my @p = split '\t', $_;
        $sample_siz{$p[0]} = $p[1]-1; }
    close SIZ; }


my %sample_grp = my %all_grp = my @agrp = ();
open SAM, "$isam" or die "Can't open $isam!";
while (<SAM>) {
    if($_ =~ /^#/) { next; }
    if(($_ =~ /sample/) or ($_ =~ /gs[em]/)) { next; }
    chomp($_); my @p = split '\t', $_;
    if($ifil) { unless(exists $incl_samp{$p[0]}) { next; } }
    
    my @q = split '\|', $p[1];
    foreach my $grp (@q) {
        $sample_grp{$p[0]}{$grp}++;
        unless(exists $all_grp{$grp}) {
            push(@agrp, $grp);
            $all_grp{$grp}++; } } }
close SAM;


my %genes = ();
if($igen) {
    open GEN, "$igen" or die "Can't open $igen!";
    while (<GEN>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $genes{$_}++; }
    close GEN; }


my $l = 0; my %idx_samp = my %incl_grp = ();
open IMAT, "$imat" or die "Can't open $imat!";

while (<IMAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    if($l == 0) {
        for(my $j=1+$iskip; $j<=$#p; $j++) {
            unless(exists $sample_grp{$p[$j]}) { next; }
            foreach my $grp (keys %{$sample_grp{$p[$j]}}) {
                $incl_grp{$grp}{$p[$j]}++; }
            $idx_samp{$j} = $p[$j]; }
        $l++;

        open OMAT, ">$omat"; print OMAT "$p[0]";
        if($ovar) { open OVAR, ">$ovar"; print OVAR "$p[0]"; }
        for(my $j=1; $j<=$iskip; $j++) {
            print OMAT "\t$p[$j]";
            if($ovar) {
                print OVAR "\t$p[$j]"; } }

        if($ogrp) { open OGRP, ">$ogrp"; print OGRP "group\tnsam\tsamples\n"; }
        foreach my $grp (@agrp) {
            unless(exists $incl_grp{$grp}) { next; }
            print OMAT "\t$grp";
            if($ovar) { print OVAR "\t$grp"; }
            if($ogrp) {
                print OGRP "$grp\t", scalar keys %{$incl_grp{$grp}};
                print OGRP "\t", join " ", keys %{$incl_grp{$grp}};
                print OGRP "\n"; } }
        print OMAT "\n";
        if($ovar) { print OVAR "\n"; }
        if($ogrp) { close OGRP; }

        next; }

    if($igen) {
        unless(exists $genes{$p[0]}) {
            next; } }

    my %gene_exp = my %gene_siz = ();
    for(my $j=1+$iskip; $j<=$#p; $j++) {
        unless(exists $idx_samp{$j}) { next; }
        foreach my $grp (keys %{$sample_grp{$idx_samp{$j}}}) {
            push(@{$gene_exp{$grp}}, $p[$j]);
            if($isiz) {
                push(@{$gene_siz{$grp}}, $sample_siz{$idx_samp{$j}}); } } }

    if($icut) {
        my ($maxm, $varm);
        if($isiz) {
            ($maxm, $varm) = get_summary($gene_exp{$agrp[0]}, $gene_siz{$agrp[0]}); }
        else {
            ($maxm, $varm) = get_summary($gene_exp{$agrp[0]}); }
        foreach my $grp (@agrp) {
            my ($m, $v);
            if($isiz) {
                ($m, $v) = get_summary($gene_exp{$grp}, $gene_exp{$grp}); }
            else {
                ($m, $v) = get_summary($gene_exp{$grp}); }
            if($m > $maxm) { $maxm = $m; } }

        if($maxm < $icut) { next; } }

    print OMAT "$p[0]";
    if($ovar) { print OVAR "$p[0]"; }
    for(my $j=1; $j<=$iskip; $j++) {
        print OMAT "\t$p[$j]"; }
    foreach my $grp (@agrp) {
        unless(exists $incl_grp{$grp}) { next; }
        my ($m, $v);
        if($isiz) {
            #print "$p[0]\t$grp\n";
            ($m, $v) = get_summary($gene_exp{$grp}, $gene_siz{$grp}); }
        else {
            ($m, $v) = get_summary($gene_exp{$grp}); }

        if($m == 0) { print OMAT "\t$m"; }
        else { printf OMAT "\t%.6f", $m; }

        if($ovar) {
            if($v == 0) { print OVAR "\t$v"; }
            else { printf OVAR "\t%.6f", $v; } } }
    print OMAT "\n";
    if($ovar) { print OVAR "\n"; } }
close IMAT; close OMAT;
if($ovar) { close OVAR; }


# Calculated mean/median of an array
sub get_summary {
    my $aref = shift;
    my $sref = shift;

    if($istat eq 'mean') {
        my $mean = my $var = 0;
        if($sref) {
            my $wsiz = my $wsum = 0;
            for(my $i=0; $i<=$#{$aref}; $i++) {
                #print "\t${$aref}[$i]\t${$sref}[$i]\n";
                $wsiz += ${$sref}[$i];
                $wsum += ${$sref}[$i]*${$aref}[$i]; }
            #print "\twsum: $wsum\twsiz: $wsiz\n"; exit;
            if($wsiz>0) {
                $mean = $wsum/$wsiz;
                $var = 0;
                for(my $i=0; $i<=$#{$aref}; $i++) {
                    $var += ${$sref}[$i]*((${$aref}[$i] - $mean)**2); }
                $var /= $wsiz; } }
        else {
            my $omean = my $n = 0;
            foreach (@$aref) {
                $n++;
                $omean = $mean;
                $mean += ($_ - $omean)/$n;
                $var += ($_ - $omean)*($_ - $mean); }
            if($n>1) { $var /= ($n-1); } }
        return ($mean, $var); }

    elsif($istat eq 'median') {
        my $median = my $mad = 0;
        if($sref) {
            die "\ncan't calc weighted median yet!\n"; }
        else {
            my @ary = sort {$a <=> $b} grep {!/^NA$/} (@{$aref});
            if(@ary % 2) { $median = $ary[@ary/2]; }
            else { $median = (($ary[(@ary/2)-1] + $ary[@ary/2]) / 2); }
            $mad = get_mad(\@ary, $median); }
        return($median, $mad); } }


sub get_mad {
    my $aref = shift;
    my $med = shift;

    my @ary = ();
    foreach (@$aref) { push(@ary, abs($_ - $med)); }
    if(@ary % 2) { return $ary[@ary/2]; }
    else { return (($ary[(@ary/2)-1] + $ary[@ary/2]) / 2); } }


__END__

=head1

Averages expression values across replicates.

=head1 USAGE

get_averaged-exp-mat.pl [--imat INPUT_MAT] [--isam sample_grp] [--istat
STATISTIC] [--icut CUTOFF] [--iskip COL_TO_SKIP] [--omat OUTPUT_MAT] [--help]

=head1 DESCRIPTION

This script takes in a gene-expression matrix and outputs a matrix with gene
values being aerages of the original values across sets of replicate samples,
specified using the sample annotations file.

=head1 ARGUMENTS

=over 12

=item C<--imat>

Input gene-expression matrix.

=item C<--isam>

Sample annotation file in tabular format with sample IDs in the first column and
original biological sample name in the second.

=item C<--icut>

(Optional) Value cut-off to print a row. Row is printed only if at least one of
the values is >= $icut.

=item C<--iskip>

(Optional) No. of columns (after the first) to skip and print as is. Default is
0.

=item C<--istat>

(Optional) 'mean' or 'median'. Default is 'mean'.

=item C<--omat>

Output matrix.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2013 May 15

=cut

