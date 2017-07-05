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

my ($help, $iobo1, $iobo2, $iobo1_filt, $ijusdef, $otab, $time, @oboref, @p);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<3);
GetOptions( 'help' => \$help,
         'iobo1=s' => \$iobo1,
     'iobo1filt=s' => \$iobo1_filt,
         'iobo2=s' => \$iobo2,
         'ijusdef' => \$ijusdef,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open TAB, ">$otab";

# Parsing OBO files
$time = runtime(); print "\n$time: Parsing the OBO files ...";

@oboref = parse_obo($iobo1);
my %obo1_id2desc = %{$oboref[0]};
my %obo1_id2syn = %{$oboref[1]};
my %obo1_id2def = %{$oboref[2]};

@oboref = parse_obo($iobo2);
my %obo2_id2desc = %{$oboref[0]};
my %obo2_id2syn = %{$oboref[1]};
my %obo2_id2def = %{$oboref[2]};

open FILT, "$iobo1_filt";
# my %words_to_avoid = ();
my %terms_tokeep = ();
while (<FILT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;
    # $words_to_avoid{$p[0]}++;
    $terms_tokeep{$p[0]}++; }
close FILT;

# Matching all pairs of terms
$time = runtime(); print "\n$time: Matching all pairs of terms ...\n";

my (%match_term, $match_dist, $all_t2, $wlength, $tmpmatch);
my (@o1d, @o2d, @tmp_dist, $wwmatch, @dist, $d, @q, @r); my $count = 0;

foreach my $term1 (keys %obo1_id2desc) {
    unless(exists $terms_tokeep{$term1}) { next; }
    $count++;
    # print "$count\t*$term1*\n\tdesc: *$obo1_id2desc{$term1}*\n\tdef: *$obo1_id2def{$term1}*\n";
    # next;

    %match_term = (); $match_dist = 0.5;

    $time = runtime();
    print "$time: $count: $term1\t$obo1_id2desc{$term1}\n";

    unless ($ijusdef) {
        @o1d = (); push(@o1d, $obo1_id2desc{$term1});
        if(exists $obo1_id2syn{$term1}) {
            push(@o1d, @{$obo1_id2syn{$term1}}) };

        # $time = runtime(); print "$time: $term1\n";
        foreach my $t1 (@o1d) {
            $t1 =~ s/[;'()]//g; $t1 =~ s/[\-\/]/ /g;
            print "\t\"$t1\"\n"; }

        foreach my $term2 (keys %obo2_id2desc) {
            @o2d = (); push(@o2d, $obo2_id2desc{$term2});
            if(exists $obo2_id2syn{$term2}) { push(@o2d, @{$obo2_id2syn{$term2}}) };
            if(exists $obo2_id2def{$term2}) { push(@o2d, $obo2_id2def{$term2}) };

            foreach (@o2d) { $_ =~ s/[;'()]//g; $_ =~ s/[\-\/]/ /g; }
            $all_t2 = join ' ', @o2d;

            @dist = ();

            foreach my $t1 (@o1d) {
                $wwmatch = 0; $wlength = 0;
                @q = split ' ', $t1;

                foreach my $w (@q) {
                    # if(exists $words_to_avoid{$w}) { next; }
                    $wlength++;
                    if(length($w) < 4) {
                        $tmpmatch = 0;
                        @r = split ' ', $all_t2;
                        foreach my $word (@r) {
                            if($word eq $w) { $tmpmatch = 1; } }
                        if($tmpmatch == 1) { $wwmatch++; } }
                    elsif($all_t2 =~ /$w/) { $wwmatch++; }
                    else {
                        $wwmatch += (1 - abs adistr($w, $all_t2)); } }

                if($wlength == 0) { push(@dist, 0); }
                # elsif($wlength <= 2) {
                #    if($wwmatch == $wlength) { push(@dist, 1); }
                #    else { push(@dist, 0); } }
                else {
                    push(@dist, $wwmatch/$wlength); } }

            $d = max_array(\@dist);
            # if($d < $match_dist) {
            if($d >= $match_dist) {
                $match_dist = $d;
                # push(@{$match_term{$d}}, $term2);
                $match_term{$d}{$term2}++; } } }
    
    # if($match_dist < 1) {
    #     if(exists $obo1_id2def{$term1}) {
    #         @q = split ' ', $obo1_id2def{$term1};

    #         foreach my $term2 (keys %obo2_id2desc) {
    #             $wwmatch = 0; $wlength = 0;

    #             foreach my $w (@q) {
    #                 $wlength++;
    #                 if($obo2_id2def{$term2} =~ /$w/) { $wwmatch++; }
    #                 else {
    #                     $wwmatch += (1 - abs adistr($w, $obo2_id2def{$term2})); } }

    #             # if($wlength == 0) { $d = 0; }
    #             # else {
    #             #     $d = ($wwmatch/$wlength); }
    #             if($wlength == 0) { push(@dist, 0); }
    #             else {
    #                 push(@dist, $wwmatch/$wlength); }

    #             $d = max_array(\@dist);
    #             if($d >= $match_dist) {
    #                 $match_dist = $d;
    #                 $match_term{$d}{$term2}++; } } } }

    if($match_dist < 0.6) { next; }

    foreach my $term2 (keys %{$match_term{$match_dist}}) {
        print TAB "$term1\t$obo1_id2desc{$term1}\t";
        print TAB "$term2\t$obo2_id2desc{$term2}\t$match_dist\n"; }
}

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

    return $min;
}

# Get minimum of array
sub max_array {
    my $aref = shift;
    my $max = shift(@{$aref});
    foreach (@{$aref}) { if($_ > $max) { $max = $_; } }

    return $max;
}


__END__

=head1 NAME

map_ontologies.pl - Map ontology 1 to matching terms in ontology 2.

=head1 USAGE

map_ontologies.pl [--iobo1 INPUT_ONTOLOGY1_OBO] [--iobo2 INPUT_ONTOLOGY2_OBO]
[--o OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes two ontologies and matches terms to ontology 1 to their best
matches in ontology 2. Supports approximate text matching. TODO: Support
substitution with appropriate adjectives (e.g. 'kidney' with 'renal') while
searching.

=head1 ARGUMENTS

=over 12

=item C<--iobo1>

Ontology 1 in OBO format.

=item C<--iobo1filt>

List of words to be filtered from terms in ontology 1.

=item C<--iobo2>

Ontology 2 in OBO format.

=item C<--ijusdef>

Use only definitions in obo2 to map terms in obo1.

=item C<--otab>

Output file: Pairs of matched terms from the two ontologies: <Term_obo1>
<Term_obo2>

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Apr 11

=cut
