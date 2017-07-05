#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
use Data::Dumper;
#use Time::SoFar qw(runtime);

# Matrix, Rows.Of.Interest, Cols.Of.Interest, SubMatrix
my ($f1, $f2, $f3, $f4) = @ARGV;
open MH, "$f1" or die "Can't open $f1!"; chomp(my @mat=<MH>); close MH;
open RH, "$f2" or die "Can't open $f2!"; chomp(my @rows=<RH>); close RH;
open CH, "$f3" or die "Can't open $f3!"; chomp(my @cols=<CH>); close CH;
open HH, ">$f4";

my %rows2incl = (); foreach (@rows) { $rows2incl{$_}++; }
my %cols2incl = (); foreach (@cols) { $cols2incl{$_}++; }

my ($num_rows2skip, $num_cols2skip) = (0, 2); # Excl. header row and first col

my @fields = split('\t', $mat[0]); my @idx_cols2incl = ();
for(my $j=($num_cols2skip+1); $j<=$#fields; $j++) {
    if(exists $cols2incl{$fields[$j]}) {
        push(@idx_cols2incl, $j);
    }
}

for(my $i=0; $i<=$#mat; $i++) {
    @fields = split('\t', $mat[$i]);

    if(($i<=$num_rows2skip) or (exists $rows2incl{$fields[0]})) {
        print HH "$fields[0]";
        
        my $j = 1;
        while($j <= $num_cols2skip) { print HH "\t$fields[$j]"; $j++; }

        foreach my $j (@idx_cols2incl) { print HH "\t$fields[$j]"; }
        print HH "\n";
    }
}

close HH;

