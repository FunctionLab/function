#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my ($help, $ifile, $isearch, $inmatch, $ieq, $inoh, $ofile);
my $iopr = 'NA'; my $icol = 'all';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 4);
GetOptions( 'help' => \$help,
         'ifile=s' => \$ifile,
       'isearch=s' => \$isearch,
       'inmatch=i' => \$inmatch,
          'icol=s' => \$icol,
          'iopr=s' => \$iopr,
             'ieq' => \$ieq,
            'inoh' => \$inoh,
         'ofile=s' => \$ofile) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

# open FH,"$ifile"; chomp(my @f=<FH>); close FH;
open HH,">$ofile";

my %ref = ();
open GH, "$isearch" or die "Can't open $isearch!";
while (<GH>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    $ref{$p[0]}++; }
close GH;

#unless($inoh) { print HH "$f[0]\n"; shift @f; }
my $nall = my $nselc = 0;

my @cols = ();

if(looks_like_number($icol)) {
    $iopr = 'NA';

    open FH, "$ifile" or die "Can't open $ifile!";
    while (<FH>) {
        chomp($_);
        if($_ =~ /^#/) { print HH "$_\n"; next; }
        if($_ =~ /^\t/) { print HH "$_\n"; next; }

        $nall++;
        my @p = split '\t', $_;

        if($p[$icol] =~ /^$/) { next; }

        if($ieq) {
            if(exists $ref{$p[$icol]}) {
                $nselc++; print HH "$_\n"; } }
        else {
            foreach my $pat (keys %ref) {
                if($p[$icol] =~ /$pat/) {
                    $nselc++; print HH "$_\n"; } } } }
    close FH; }
else {
    if($icol =~ /:/) {
        my @temp_cols = split '_', $icol;

        foreach my $entry (@temp_cols) {
            if($entry =~ /:/) {
                my @p = split ':', $entry;

                for(my $k=$p[0]; $k<=$p[1]; $k++) {
                    push(@cols, $k); } }
            else { push(@cols, $entry); } } }
    else { @cols = split '_', $icol; }
    print "@cols\n";

    open FH, "$ifile" or die "Can't open $ifile!";
    while (<FH>) {
        chomp($_);
        if($_ =~ /^#/) { print HH "$_\n"; next; }
        if($_ =~ /^\t/) { print HH "$_\n"; next; }
        $nall++;
        my @p = split '\t', $_;

        if($icol eq 'all') {
            @cols = ();
            for(my $j=0; $j<=$#p; $j++) {
                push(@cols, $j); } }

        if($iopr eq 'AND') {
            my $par = 1;

            if($ieq) {
                foreach my $c (@cols) {
                    unless(exists $ref{$p[$c]}) { $par = 0; last; } } }
            else {
                COL: foreach my $c (@cols) {
                    foreach my $pat (keys %ref) {
                        if($p[$c] =~ /$pat/) { next COL; } }
                    $par = 0; last COL; } }

            if($par == 1) { $nselc++; print HH "$_\n"; } }

        elsif($iopr eq 'OR') {
            if($ieq) {
                foreach my $c (@cols) {
                    if(exists $ref{$p[$c]}) {
                        $nselc++; print HH "$_\n"; last; } } }
            else {
                foreach my $c (@cols) {
                    foreach my $pat (keys %ref) {
                        if($p[$c] =~ /$pat/) { 
                            $nselc++; print HH "$_\n"; last; } } } } }

        elsif($inmatch) {
            my $tnm = 0;

            if($ieq) {
                foreach my $c (@cols) {
                    if(exists $ref{$p[$c]}) { $tnm++; } } }
            else {
                COL: foreach my $c (@cols) {
                    foreach my $pat (keys %ref) {
                        if($p[$c] =~ /$pat/) { $tnm++; next COL; } } } }

            if($tnm >= $inmatch) { $nselc++; print HH "$_\n"; } } }
    close FH; }


print "\nTot. no. lines in file: $nall\nNo. lines selected: $nselc\n\n";

close HH;


__END__

=head1

Greps input file using the search patterns listed in a file and returns relevant
rows. Does more than 'grep -f' because specific columns and logical operators
can be specified.

=head1 USAGE

grep_minusf.pl [--ifile INPUT_FILE] [--isearch INPUT_PATTERNS] [--icol
COLUMNS_TO_INCL] [--iopr LOGICAL_OPERATOR] [--ofile OUTPUT_FILE] [--help]

=head1 DESCRIPTION

This script takes in a file and a list of search patterns (in a file) and
returns rows in the file that contain any one of the patterns. Search can be
customized by specifying specific columns of the file to search for and logical
operaqtors to combine the search across columns.

=head1 ARGUMENTS

=over 12

=item C<--ifile>

Input file in tab-delimited format.

=item C<--isearch>

List of strings to search for.

=item C<--icol>

Column(s) to include in search. Index begins at 0; Multiple indices can be
provided by stringing them together using underscores. Column ranges can be
provided using ":". For e.g., 3_5:7 means columns 3,5,6,7.

=item C<--iopr>

(Optional) NA, AND or OR. Relevant when choosing more than one column as the
basis of filtering. Signifies whether to look for presence of reference in one
of or all columns. Default 'NA'.

=item C<--ieq>

(Optional) If provided, "column_entry 'equals' pattern" will be checked instead
of "column_entry 'contains' pattern".

=item C<--inoh>

(Optional) If provided, the first line of the input file is assumed to NOT
contain the header/column-names and is used in search.

=item C<--ofile>

Output file containing one of the seach strings satisfying the
column/logical-operator specifications.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

