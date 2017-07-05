#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

my($help, $in_ga, $out_ga); my @ec = (); my $ns = '';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(  'help' => \$help,
            'iga=s' => \$in_ga,
          'ec=s{,}' => \@ec,
             'ns=s' => \$ns,
            'oga=s' => \$out_ga) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if (($help) or (($#ec == -1) and ($ns eq '')));


open GA, "$in_ga" or die "Can't open file: $in_ga!"; chomp(my @ga=<GA>); close GA;
open FH, ">$out_ga";

foreach (@ga) {
    if($_ =~ /^!/) { print FH "$_\n"; } }

my (@p, $par, %ecs); my $tot_ann = my $num_sel_ann = 0;

if($#ec >= 0) {
    %ecs = (); foreach(@ec) { $ecs{$_}++; } }

foreach (@ga) {
    unless($_ =~ /^!/) {
        $tot_ann++;
        @p = split '\t', $_;

        $par = 0;
        if(($#ec >= 0) and ($ns ne '')) {
            if((exists $ecs{$p[6]}) and ($p[8] eq $ns)) { $par = 1; }
        }
        else {
            if(($#ec >= 0) and (exists $ecs{$p[6]})) { $par = 1; }
            if(($ns ne '') and ($p[8] eq $ns)) { $par =1; }
        }
        if ($par == 1) { print FH "$_\n"; $num_sel_ann++; }
    }
}

print "\nTot. ann in $in_ga: $tot_ann\nNo. of filtered ann: $num_sel_ann\n\n";

close FH;

__END__

=head1 NAME

filter_gene_association_file.pl
	- trims GO gene association file to specified evidences or name_space

=head1 USAGE

./filter_gene_association_file_by_evidenccode.pl [--iga INPUT_GAFILE] [--ec EVIDENCE_CODES] [--ns NAME_SPACE] [--oga OUTPUT_GA] [--help]

=head1 DESCRIPTION

This script takes a typical gene_association file and filters it to include only
annotations based on a use-specified set of evidence codes or names_space.

=head1 ARGUMENTS

=over 12

=item C<--iga>

Annotation file - gene_association file - obtained from GO

=item C<--ec>

Evidence codes, e.g. IPI IEP ISS ... These have to be the 3-letter codes (separated by spaces) specified in http://www.geneontology.org/GO.evidence.shtml

=item C<--ns>

Name_space: One of P, F or C

=item C<--oga>

Output annotation file in the same format, only containing annotaions based on
chosen ECs or Name_space

=item C<--help>

prints this documentation

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2011 October 12

=cut

