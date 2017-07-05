#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
#use String::Approx qw(adist adistr);
use Time::SoFar qw(runtime);

sub parse_obo;

my ($help, $iobo, $ifilt, $otab, $time, @p);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<3);
GetOptions( 'help' => \$help,
          'iobo=s' => \$iobo,
         'ifilt=s' => \$ifilt,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


# Parsing OBO files
$time = runtime(); print "\n$time: Parsing the OBO files ...";

my @oboref = parse_obo($iobo);
my %obo_id2desc = %{$oboref[0]};
my %obo_id2syn = %{$oboref[1]};
my %obo_id2def = %{$oboref[2]};


my %filtt = ();
if($ifilt) {
# Parsing list of filtred terms to retain
    $time = runtime(); print "\n$time: Parsing filtered terms ...";

    open FIL, "$ifilt" or die "Can't open $ifilt!";
    while (<FIL>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $filtt{$_}++; }
    close FIL; }


# Printing list of all terms
$time = runtime(); print "\n$time: Printing list of all terms ...\n";

open TAB, ">$otab";

foreach my $term (keys %obo_id2desc) {
    if($ifilt) { unless(exists $filtt{$term}) { next; } }

    print TAB "$term\t$obo_id2desc{$term}\tname\n";
    if(exists $obo_id2def{$term}) {
        print TAB "$term\t$obo_id2def{$term}\tdef\n"; }

    if(exists $obo_id2syn{$term}) {
        foreach my $syn (@{$obo_id2syn{$term}}) {
            print TAB "$term\t$syn\tsyn\n"; } } }


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
            #shift(@p); splice(@p, -2);
            ($syn = $line) =~ s/^synonym: \"([^\"]*)\" [A-Z].*$/$1/g;
            #$syn = lc(join ' ', @p);
            $syn =~ s/\"//g;
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


__END__

=head1

Extract list of terms and descriptions from ontology.

=head1 USAGE

get_list-from-obo.pl [--iobo ONTOLOGY_OBO] [--otab OUTPUT_LIST] [--help]

=head1 DESCRIPTION

This script takes an ontology and extracts the list of ids and terms including
definitions and synonymns. TODO: Support substitution with appropriate
adjectives (e.g. 'kidney' with 'renal') while searching.

=head1 ARGUMENTS

=over 12

=item C<--iobo>

Ontology in OBO format.

=item C<--otab>

List of terms and available descriptive tags from the ontology.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Apr 11

=cut
