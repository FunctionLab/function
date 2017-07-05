#!/usr/bin/perl
use strict;
use warnings;
use Graph;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my($help, $in_ann, $in_obo, $in_map, $out_ann, $time); my $in_slim = '';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<4);
GetOptions( 'help' => \$help,
          'iann=s' => \$in_ann,
          'iobo=s' => \$in_obo,
         'islim=s' => \$in_slim,
          'imap=s' => \$in_map,
          'oann=s' => \$out_ann) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open ANN, "$in_ann" or die "Can't open $in_ann!"; chomp(my @ann=<ANN>); close ANN;
open OBO, "$in_obo" or die "Can't open $in_obo!"; chomp(my @obo=<OBO>); close OBO;
open MAP, "$in_map" or die "Can't open $in_map!"; chomp(my @map=<MAP>); close MAP;
open PRO, ">$out_ann";

# Parson OBO file
$time = runtime(); print "\n$time: Parsing the OBO file ...";

my $inside = 0; my (@p, $term, $name, $syn, $rel, $parent);
my %id_2desc = (); my %alt2std_id = (); # my %desc_2id = ();
# my %syn_desc_2id = (); my %syn_count = ();

my %is_a = (); my %part_of = (); #%regulates = ();
my $dag = Graph->new;

LINE: foreach my $line (@obo) {
    if($line =~ /^#/) { next; }
    @p = split(' ', $line);

    if ($line =~ /^$/) { next LINE; }
    elsif ($p[0] eq '[Term]') { $inside = 1; }
    elsif ($inside and ($p[0] eq 'id:')) { $term = $p[1]; }
    elsif ($inside and ($p[0] eq 'name:')) {
        shift(@p); $name = lc(join ' ', @p);
        $id_2desc{$term} = $name;
    }
    elsif ($inside and ($p[0] eq 'alt_id:')) {
        # push(@{$alt2std_id{$p[1]}}, $term);
        $alt2std_id{$p[1]} = $term;
    }
    # elsif ($inside and ($p[0] eq 'synonym:')) {
    #    shift(@p); splice(@p, -2);
    #    $syn = lc(join ' ', @p); $syn =~ s/\"//g;
    #    $syn_count{$syn}++; $syn_desc_2id{$syn} = $term;
    # }
    elsif ($inside and ($p[0] eq 'is_a:')) {
        $parent = splice(@p, 0, 2);
        $is_a{join '__', ($term, $parent)}++;
        $dag->add_edge($term, $parent);
    }
    elsif ($inside and ($p[0] eq 'relationship:')) {
        $rel = $p[1]; $parent = splice(@p, 0, 3);
        
        if($rel eq 'has_part') { next LINE; }
        elsif($rel eq 'part_of') {
            $part_of{join '__', ($term, $parent)}++;
            $dag->add_edge($term, $parent);
        }
    }
    elsif ($inside and ($p[0] eq 'is_obsolete:')) {
        delete $id_2desc{$term};
    }
}

# Remove synonymns that match more than one term
# foreach my $syn (keys %syn_desc_2id) {
#    if($syn_count{$syn} > 1) { delete $syn_desc_2id{$syn}; }
# }

my %slim_ids = ();
unless($in_slim eq '') {
    open SLIM, "$in_slim" or die "Can't open $in_slim!";
    chomp(my @slim=<SLIM>); close SLIM;

    foreach (@slim) {
        @p = split '\t', $_;
        unless(exists $id_2desc{$p[0]}) { next; }
        $slim_ids{$p[0]} = $p[1];
    }
}

# Mapping given gene ids to new ids
$time = runtime(); print "\n$time: Mapping given gene ids to new ids ...";

my %id_map = ();
foreach (@map) {
    @p = split '\t', $_;
    $id_map{$p[0]} = $p[1];
}

# Parsing annotation file
$time = runtime(); print "\n$time: Parsing gene annotations ...";

my %term_genes = ();
foreach (@ann) {
    @p = split '\t', $_;
    unless(exists $id_map{$p[0]}) { next; }

    if(exists $id_2desc{$p[1]}) {
        $term_genes{$p[1]}{$id_map{$p[0]}}++;
    }
    elsif(exists $alt2std_id{$p[1]}) {
        $term_genes{$alt2std_id{$p[1]}}{$id_map{$p[0]}}++;
    }
}

# Topological sorting and propagation
$time = runtime(); print "\n$time: Topologically sorting all the terms ...";
my @ts_terms = $dag->topological_sort;

# Propagating gene annotations
$time = runtime(); print "\n$time: Propagating annotations ...";

# my (@path, @term_slim, $slim_tag);
my $num_terms = 0; my %final_genes = ();
my %prop_term_genes = (); my ($num_genes, @parents);
my (@ancestors, %anc_idx, @term_slim, $slim_tag);
foreach my $term (@ts_terms) {
    if(exists $term_genes{$term}) {
        foreach my $gene (keys %{$term_genes{$term}}) {
            $prop_term_genes{$term}{$gene}++;
        }
    }

    unless(exists $prop_term_genes{$term}) { next; }
    print PRO "$term\t$id_2desc{$term}"; $num_terms++;

    # Assigning slim terms if this option is selected
    unless($in_slim eq '') {
        @ancestors = all_successors($term);
        foreach (@ancestors) { $anc_idx{$_}++; }
        @term_slim = ();
        foreach my $sterm (keys %slim_ids) {
            # @path = $dag->SP_Dijkstra($term, $sterm);
            # unless(scalar(@path) == 0) {
            #     push(@term_slim, $sterm);
            # }
            if(exists $anc_idx{$sterm}) {
                push(@term_slim, $sterm);
            }
        }
        $slim_tag = join '|', @term_slim;
        # if(scalar(@term_slim) == 0) { $slim_tag = 'MP:0005395'; }
        print PRO " $slim_tag";
    }

    $num_genes = scalar(keys %{$prop_term_genes{$term}});
    print PRO " ($num_genes)";

    foreach my $gene (keys %{$prop_term_genes{$term}}) {
        $final_genes{$gene}++;
        print PRO "\t$gene";

        @parents = $dag->successors($term);
        foreach my $pt (@parents) {
            $prop_term_genes{$pt}{$gene}++;
        }
    }
    print PRO "\n";
}

$num_genes = scalar(keys %final_genes);

$time = runtime(); print "\nDONE $time\n\n";
print "No. Terms: $num_terms\nNo. Genes: $num_genes\n\n";

close PRO;

__END__

=head1

Propagates annotations using the true-path-rule

=head1 USAGE

./propagate_annotations_upaDAG.pl [--iann DIRECT_ANNOTATIONS] [--iobo ONTOLOGY_OBO]
[--imap GENE_IDMAP] [--oann PROPAGATED_ANNOTATIONS] [--help]

=head1 DESCRIPTION

This script takes (i) direct annotations of genes to terms in an ontology, &
(ii) ontology of terms and their relationships, to gives annotations by
propagating annotations up the ontology based on 'is-a' and 'part-of'
relationships between terms.

=head1 ARGUMENTS

=over 12

=item C<--iann>

File containing direct annotations of genes to terms: <GeneID> <TermID>
tab-delimited; one/line.

=item C<--iobo>

Ontology file in OBO format

=item C<--islim>

(Optional) Slim terms in the ontology: <TermID> <Desc> tab-delimited; one/line.
If this file is provided, then, (i) an additional output file containing mapping of all
terms to their slim terms is created, and (ii) the slim term corresponding to a
given term is printed alongside the term description in the GMT file.

=item C<--imap>

Gene ID mapping file: <GeneID_old> <GeneID_new> tab-delimited; one/line

=item C<--oann>

Output file containing gene annotations propagated to terms in the ontology.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Feb 25

=cut


