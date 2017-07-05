#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
#use lib '/Genomics/ogtr04/arjunk/src/perl/downloads/';
use String::Approx qw(adist adistr);
use Time::SoFar qw(runtime);

sub min_array;
sub max_array;

#my ($help, $ilist1, $ilist2, $itype, $ionly1to2, $otab1, $otab2, $time, @oboref, @p);
my ($help, $ilist1, $ilist2, $itype, $inosplit, $ionly1to2, $otab, $time, @oboref, @p);
my $icut = 0.9; my $iallp = 1; my $icomb = 'max';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<3);
GetOptions( 'help' => \$help,
        'ilist1=s' => \$ilist1,
        'ilist2=s' => \$ilist2,
         'itype=s' => \$itype,
       'ionly1to2' => \$ionly1to2,
        'inosplit' => \$inosplit,
         'icomb=s' => \$icomb,
          'icut=f' => \$icut,
         'iallp=s' => \$iallp,
         #'otab1=s' => \$otab1,
         #'otab2=s' => \$otab2,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

# open TAB1, ">$otab1";
# open TAB2, ">$otab2";
open TAB, ">$otab";

# Parsing LIST and OBO files
$time = runtime(); print "\n$time: Parsing lists ...";

open LIST, "$ilist1";
my %list1_terms = my %list1_type = ();
while (<LIST>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    if($itype) { if($p[2] ne $itype) { next; } }
    $list1_terms{$p[1]} = $p[0];
    $list1_type{$p[1]} = $p[2]; }
close LIST;

open LIST, "$ilist2";
my %list2_terms = my %list2_type = ();
while (<LIST>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    if($itype) { if($p[2] ne $itype) { next; } }
    $list2_terms{$p[1]} = $p[0];
    $list2_type{$p[1]} = $p[2]; }
close LIST;

# Matching all pairs of terms
$time = runtime(); print "\n$time: Matching all pairs of terms ...\n";

my ($wlength, $tmpmatch);
my (@od, @tmp_dist, $wwmatch, @dist, $d, @q, @r); my $count = 0;
my %pair_desc = my %pair_sim = ();
my %adistr = ();

T1: foreach my $ot1 (keys %list1_terms) {
    $count++;
    my %match_term = (); my $max_sim = 0.5;

    unless($count % 100) {
        $time = runtime();
        print "\t$time: $count - $ot1\n"; }
    #$time = runtime();
    #print "\t$time: $count – $ot1\n";

    my $t1 = lc($ot1);
    $t1 =~ s/[.;'()\[\],]//g; $t1 =~ s/[\-\/_]/ /g;
    $t1 =~ s/\" \[.*$//g; $t1 =~ s/\.$//g; $t1 =~ s/  */ /g;

    T2: foreach my $ot2 (keys %list2_terms) {
        my $t2 = lc($ot2);
        $t2 =~ s/[.;'()\[\],]//g; $t2 =~ s/[\-\/_]/ /g;
        $t2 =~ s/\" \[.*$//g; $t2 =~ s/\.$//g; $t2 =~ s/  */ /g;
        #if($t2 =~ /gsm901657 /) { print "gsm901657\t$t2\n"; }

        my $d1 = my $d2 = $wwmatch = $wlength = 0;

        if($ionly1to2 and $inosplit) {
            if(length($t1) <= 5) {
                $tmpmatch = 0;
                @r = split ' ', $t2;
                foreach my $word (@r) {
                    #if($t2 =~ /gsm901657 /) { print "\t\t$w --> $word\n"; }
                    if(($word eq $t1) or ($word.'s' eq $t1) or ($word eq $t1.'s')) { $tmpmatch = 1; } }
                if($tmpmatch == 1) { $d1 = 1; }
                #if($t2 =~ /gsm901657 /) { print "\t$t1\t$wlength\t$w --> $tmpmatch --> $wwmatch\n"; }
            }
            #elsif($t2 =~ /$t1/) { $d1 = 1; }
            else {
                @q = split ' ', $t1;
                my $shortwnum = my $shortwmat = 0;
                foreach my $w (@q) {
                    $tmpmatch = 0; if(length($w) <= 5) { $shortwnum++; }
                    @r = split ' ', $t2;
                    foreach my $word (@r) {
                        if(($word eq $w) or ($word.'s' eq $w) or ($word eq $w.'s')) {
                            $tmpmatch = 1; } }
                    if($tmpmatch == 1) {
                        $wwmatch++;
                        if(length($w) <= 5) { $shortwmat++; } } }

                if(($wwmatch > 0) and ($shortwmat eq $shortwnum)) {
                    my $wpair = join '__', sort($t1, $t2);
                    unless(exists $adistr{$wpair}) {
                        $adistr{$wpair} = abs adistr($t1, $t2); }
                    $d1 += (1 - $adistr{$wpair}); } } }

        unless($inosplit) {
            @q = split ' ', $t1;
            foreach my $w (@q) {
                $wlength++;
                #if($t2 =~ /gsm901657 /) { print "\t$t1\t$wlength\t$w\n"; }
                if(length($w) <= 5) {
                    $tmpmatch = 0;
                    @r = split ' ', $t2;
                    foreach my $word (@r) {
                        #if($t2 =~ /gsm901657 /) { print "\t\t$w --> $word\n"; }
                        if(($word eq $w) or ($word.'s' eq $w) or ($word eq $w.'s')) { $tmpmatch = 1; } }
                    if($tmpmatch == 1) { $wwmatch++; }
                    #if($t2 =~ /gsm901657 /) { print "\t$t1\t$wlength\t$w --> $tmpmatch --> $wwmatch\n"; }
                }
                elsif($t2 =~ /$w/) { $wwmatch++; }
                else {
                    my $wpair = join '__', sort($w, $t2);
                    unless(exists $adistr{$wpair}) {
                        $adistr{$wpair} = abs adistr($w, $t2); }
                    $wwmatch += (1 - $adistr{$wpair}); } }

            if($wlength > 0) {
                $d1 = ($wwmatch/$wlength); } }

        if($ionly1to2) { $d = $d1; }
        else {
            @q = split ' ', $t2; # print "\n$t2\n";
            $wwmatch = $wlength = 0;

            foreach my $w (@q) {
                $wlength++; # print "\t$w\n";
                if(length($w) <= 5) {
                    $tmpmatch = 0;
                    @r = split ' ', $t1;
                    foreach my $word (@r) {
                        if($word eq $w) { $tmpmatch = 1; } }
                    if($tmpmatch == 1) { $wwmatch++; } }
                elsif($t1 =~ /$w/) { $wwmatch++; }
                else {
                    my $wpair = join '__', sort($w, $t1);
                    unless(exists $adistr{$wpair}) {
                        $adistr{$wpair} = abs adistr($w, $t1); }
                    $wwmatch += (1 - $adistr{$wpair}); } }

            if($wlength > 0) {
                $d2 = ($wwmatch/$wlength); }

            $d = $d1;
            if($icomb eq 'max') {
                if($d2 > $d1) { $d = $d2; } }
            elsif($icomb eq 'min') {
                if($d2 < $d1) { $d = $d2; } } }

        if($d < $icut) { next T2; }

        # print TAB1 "$list1_terms{$ot1}\t$ot1\t$list1_type{$ot1}\t";
        # print TAB1 "$list2_terms{$ot2}\t$ot2\t$list2_type{$ot2}\t";
        # if($ionly1to2) {
        #     printf TAB1 "%.3f\n", $d; }
        # else {
        #     printf TAB1 "%.3f\t%.3f\t%.3f\n", $d1, $d2, $d; }


        my $pair = join '___', ($list1_terms{$ot1}, $list2_terms{$ot2});

        if($iallp) {
            if(exists $pair_sim{$pair}) {
                if($d > $pair_sim{$pair}) {
                    $pair_desc{$pair} = join '___', ($ot1, $ot2);
                    $pair_sim{$pair} = $d; } }
            else {
                $pair_desc{$pair} = join '___', ($ot1, $ot2);
                $pair_sim{$pair} = $d; } }
        else {
            if($d >= $max_sim) {
                $max_sim = $d;
                $match_term{$d}{$ot2}++; } } }

    # if($max_sim < 0.5) { next; }

    unless($iallp) {
        foreach my $ot2 (keys %{$match_term{$max_sim}}) {
            print TAB "$list1_terms{$ot1}\t$ot1\t$list1_type{$ot1}\t";
            print TAB "$list2_terms{$ot2}\t$ot2\t$list2_type{$ot2}\t";
            printf TAB "%.3f\n", $max_sim; } } }

#close TAB1;

if($iallp) {
# Printing best matches
    foreach my $pair (sort keys %pair_sim) {
        my ($id1, $id2) = split '___', $pair;
        my ($ot1, $ot2) = split '___', $pair_desc{$pair};
        print TAB "$id1\t$ot1\t$list1_type{$ot1}\t";
        print TAB "$id2\t$ot2\t$list2_type{$ot2}\t";
        printf TAB "%.3f\n", $pair_sim{$pair}; } }

# close TAB2;
close TAB;

$time = runtime(); print "\n$time: DONE\n\n";


# Subroutines
# ===========
# Get minimum of array
sub min_array {
    my $aref = shift;
    my $min = shift(@{$aref});
    foreach (@{$aref}) { if($_ < $min) { $min = $_; } }

    return $min; }

# Get minimum of array
sub max_array {
    my $aref = shift;
    my $max = shift(@{$aref});
    foreach (@{$aref}) { if($_ > $max) { $max = $_; } }

    return $max; }


__END__

=head1

Map 2 lists of terms to eachother.

=head1 USAGE

map_list-to-list.pl [--ilist11 INPUT_LIST1] [--ilist2 INPUT_LIST2] [--otab1
OUTPUT_TABLE1] [--otab2 OUTPUT_TABLE2] [--help]

=head1 DESCRIPTION

This script takes two lists and matches terms in list 1 to their best
matches in list 2. Supports approximate text matching. TODO: Support
substitution with appropriate adjectives (e.g. 'kidney' with 'renal') while
searching.

=head1 ARGUMENTS

=over 12

=item C<--ilist1>

First list of terms/phrases to map.

=item C<--ilist2>

Second list of terms/phrases to map.

=item C<--ionly1to2>

(Optional) Map terms only from first list to second (not second onto first). This way, for
e.g., term1 from list1 that is a substring of term2 of list2 will get a perfect
score, but term2 being a substring of term1 will not receive any score. Useful
for matching a list of names & synonyms to a list of descriptions & definitions.

=item C<--inosplit>

(Optional) Applicable only when ionly1to2 is set. Do not split term1 of list1 into
individual words, but retain and match them as phrases to terms in list2.

=item C<--icut>

(Optional) Cut-off score. Default is 0.9.

=item C<--iallp>

(Optional) Output all pairs of matches (1) or just the best matches in list2 per
term in list 1. Default 1.

=item C<--otab>

Output file: Pairs of matched terms with score > icut: <List1_Term> <List2_Term> <Score>

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Apr 11

=cut
