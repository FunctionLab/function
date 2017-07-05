#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use String::Approx qw(adist adistr);
use Time::SoFar qw(runtime);

sub parse_obo;
sub min_array;
sub max_array;

my ($help, $iobo, $ilist, $otab, $time, @oboref, @p);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<3);
GetOptions( 'help' => \$help,
          'iobo=s' => \$iobo,
         'ilist=s' => \$ilist,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open TAB, ">$otab";

# Parsing LIST and OBO files
$time = runtime(); print "\n$time: Parsing the LIST and OBO files ...";

@oboref = parse_obo($iobo);
my %obo_id2desc = %{$oboref[0]};
my %obo_id2syn = %{$oboref[1]};
my %obo_id2def = %{$oboref[2]};

open LIST, "$ilist";
my %list_terms = ();
while (<LIST>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    $list_terms{$p[0]}++; }
close LIST;

# Matching all pairs of terms
$time = runtime(); print "\n$time: Matching all pairs of terms ...\n";

my (%match_term, $match_dist, $allt, $wlength, $tmpmatch);
my ($t1, @od, @tmp_dist, $wwmatch, @dist, $d, @q, @r); my $count = 0;

foreach my $ot1 (keys %list_terms) {
    $count++;
    %match_term = (); $match_dist = 0.5;

    $time = runtime();
    print "\t$time: $count - $ot1";

    $t1 = lc($ot1);
    $t1 =~ s/[.;'()]//g; $t1 =~ s/[\-\/]/ /g;
    $t1 =~ s/\" \[.*$//g; $t1 =~ s/\.$//g;
    print " > $t1\n";

    foreach my $t2 (keys %obo_id2desc) {
        @od = (); push(@od, $obo_id2desc{$t2});
        if(exists $obo_id2syn{$t2}) { push(@od, @{$obo_id2syn{$t2}}) };
        if(exists $obo_id2def{$t2}) { push(@od, $obo_id2def{$t2}) };

        foreach (@od) { $_ =~ s/[.;'()]//g; $_ =~ s/[\-\/]/ /g; }
        $allt = join ' ', @od;

        $wwmatch = 0; $wlength = 0;
        @q = split ' ', $t1;

        foreach my $w (@q) {
            $wlength++;
            if(length($w) < 5) {
                $tmpmatch = 0;
                @r = split ' ', $allt;
                foreach my $word (@r) {
                    if($word eq $w) { $tmpmatch = 1; } }
                if($tmpmatch == 1) { $wwmatch++; } }
            elsif($allt =~ /$w/) { $wwmatch++; }
            else {
                $wwmatch += (1 - abs adistr($w, $allt)); } }

        if($wlength == 0) { $d = 0; }
        else {
            $d = ($wwmatch/$wlength); }

        if($d >= $match_dist) {
            $match_dist = $d;
            $match_term{$d}{$t2}++; } }

    if($match_dist < 0.5) { next; }

    foreach my $t2 (keys %{$match_term{$match_dist}}) {
        print TAB "$ot1\t$t2\t$obo_id2desc{$t2}\t$match_dist\n"; } }

close TAB;
$time = runtime(); print "\n$time: DONE\n\n";

# Subroutines
# ===========
# Parse OBO and return ids and descriptions
sub parse_obo {
    my $obofile = shift;
    open OBO, "$obofile" or die "Can't open $obofile!";
        chomp(my @obo=<OBO>); close OBO;

    my $in = 0; my (@p, $term, $name, $def, $syn);
    my %id2desc = (); my %id2syn = (); my %id2def = ();
    my %syn_count = (); my @temp_array;

    LINE: foreach my $line (@obo) {
        if($line =~ /^#/) { next LINE; }
        @p = split(' ', $line);

        if ($line =~ /^$/) { next LINE; }
        elsif ($p[0] eq '[Term]') { $in = 1; }
        elsif ($p[0] eq '[Typedef]') { $in = 0; }
        elsif ($in and ($p[0] eq 'id:')) { $term = $p[1]; }
        elsif ($in and ($p[0] eq 'name:')) {
            shift(@p); $name = lc(join ' ', @p);
            $id2desc{$term} = $name;
        }
        elsif ($in and ($p[0] eq 'namespace:')) {
            if(($term =~ /^GO:/) and ($line !~ /biological_process/)) {
                delete $id2desc{$term}; $in = 0; next LINE;
            }
        }
        elsif ($in and ($p[0] eq 'def:')) {
            ($def = lc($line)) =~ s/^def: \"//g;
            $def =~ s/\" \[.*$//g; $def =~ s/\.$//g;
            $def =~ s/['()]//g; $def =~ s/[\-\/]/ /g;
            $id2def{$term} = $def;
        }
        elsif ($in and ($p[0] eq 'synonym:')) {
            if(($term =~ /^GO:/) and (($line !~ /EXACT/) and ($line !~ /ADDED/))) { next LINE; }
            if(($term=~ /^BTO/) and (($line =~ /RELATED GE/) or ($line =~ /RELATED SCI/))) { next LINE; }
            shift(@p); splice(@p, -2);
            $syn = lc(join ' ', @p); $syn =~ s/\"//g;
            if(length($syn) < 5) { next LINE; }
            $syn_count{$syn}++;
            push(@{$id2syn{$term}}, $syn);
        }
        elsif ($in and ($p[0] eq 'is_obsolete:')) {
            delete $id2desc{$term};
        }
    }

    # Remove synonymns that match more than one term
    foreach my $term (keys %id2syn) {
        @temp_array =  ();
        foreach my $syn (@{$id2syn{$term}}) {
            if($syn_count{$syn} > 1) { next; }
            push(@temp_array, $syn);
        }
        if(@temp_array == 0) { delete $id2syn{$term}; }
        else { @{$id2syn{$term}} = @temp_array; }
    }

    return (\%id2desc, \%id2syn, \%id2def);
}

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

    return $max;
}


__END__

=head1

Map list of terms/phrases to matching terms in ontology.

=head1 USAGE

map_list_to_ontology.pl [--iobo ONTOLOGY_OBO] [--ilist LIST-OF-TERMS] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes a list of terms/phrases and maps them to terms in an ontology.
Supports approximate text matching. TODO: Support substitution with appropriate
adjectives (e.g. 'kidney' with 'renal') while searching.

=head1 ARGUMENTS

=over 12

=item C<--iobo>

Ontology in OBO format.

=item C<--ilist>

List of terms/phrases to map to terms in the ontology.

=item C<--ijusdef>

Use only definitions in obo to map terms in obo1. [CURRENTLY NOT IN USE.]

=item C<--otab>

Output file: Pairs of matched terms from the two ontologies:
<List_Term> <Ontology_Term> ...

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Apr 11

=cut
