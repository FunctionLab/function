#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util 'shuffle';

my ($help, $in_dat, $in_frac, $in_num, $out_test, $ofile1, $ofile2);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'idat=s' => \$in_dat,
         'ifrac=s' => \$in_frac,
          'inum=s' => \$in_num,
           'otest' => \$out_test) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open FH, "$in_dat" or die "Can't open $in_dat!";
    chomp(my @dat=<FH>); close FH;

($ofile1 = $in_dat) =~ s/\.([a-z]*)$/\.sub\.$1/g;
if($out_test) {
    ($ofile1 = $in_dat) =~ s/\.([a-z]*)$/\.train\.$1/g;
    ($ofile2 = $ofile1) =~ s/\.train\./\.test\./g; open TE, ">$ofile2";
}
open TR, ">$ofile1";

if($in_frac) { $in_num = int(scalar(@dat)*$in_frac + 0.5); }

my @rand_idx = ();
for(my $i=0; $i<=$#dat; $i++) { push(@rand_idx, $i); }
@rand_idx = shuffle(@rand_idx);

my %sub_idx = ();
foreach (@rand_idx[0..($in_num-1)]) { $sub_idx{$_}++; }

for(my $i=0; $i<=$#dat; $i++) {
    if(exists $sub_idx{$i}) {
        print TR "$dat[$i]\n"; }
    elsif($out_test) {
        print TE "$dat[$i]\n"; }
}

close TR; close TE;

__END__

=head1

Randomly split data into two.

=head1 USAGE

./rand_split_data.pl[--idat INPUT_DATA] [--ifrac FRACTION] [--help]

=head1 DESCRIPTION

This script takes a data file & fraction from the user and produces two files,
with the fist one containing 'fraction' of the data file and the other
containing the rest of the data.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input data file

=item C<--ifrac>

Fraction of data to keep in subset (e.g. train).

=item C<--inum>

No. of lines of data to keep in subset (e.g. train).

=item C<--otest>

Prints rest of the data in another file (e.g. test).

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Apr 24

=cut

