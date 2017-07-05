#!/usr/bin/perl
use strict;
use warnings;
use Graph;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);


my ($help, $iann, $iobo, $imap, $ogmt, $time, $islim, $iregp, $inoprop, $igsfilt, $isemsim);
my ($ilgs1, $ilgs2); my $icut = 0.8;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<2);
GetOptions( 'help' => \$help,
          'iann=s' => \$iann,
          'iobo=s' => \$iobo,
         'islim=s' => \$islim,
          'imap=s' => \$imap,
           'iregp' => \$iregp,
         'inoprop' => \$inoprop,
       'igsfilt=s' => \$igsfilt,
         'isemsim' => \$isemsim,
         'ilgs1=s' => \$ilgs1,
         'ilgs2=s' => \$ilgs2,
          'icut=f' => \$icut,
          'ogmt=s' => \$ogmt) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


# Parse gs-terms to filter
my %filt_gs = ();
if($igsfilt) {
    open FGS, "$igsfilt" or die "Can't open $igsfilt!";
    while (<FGS>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $filt_gs{$p[0]}++; }
    close FGS; }


# Parson OBO file
$time = runtime(); print "\n$time: Parsing the OBO file ...";
my ($obog, $id_2desc, $alt2std_id);
if($igsfilt) {
    ($obog, $id_2desc, $alt2std_id) = parse_obo($iobo, \%filt_gs); }
else {
    ($obog, $id_2desc, $alt2std_id) = parse_obo($iobo); }


# Parsing SLIM file
my %slim_ids = ();
if($islim) {
    open SLIM, "$islim" or die "Can't open $islim!";
    while (<SLIM>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        unless(exists $id_2desc->{$p[0]}) { next; }
        $slim_ids{$p[0]} = $p[1]; }
    close SLIM; }


# Topological sorting
$time = runtime(); print "\n$time: Top.sorting and printing all the terms ...";
my @ts_terms = $obog->topological_sort; # print "\n\troot: $ts_terms[$#ts_terms]\t$id_2desc->{$ts_terms[$#ts_terms]}";
my %term_level = (); my $l = 0; foreach (@ts_terms) { $term_level{$_} = $l; $l++; }
my %term_anc = my %term_nanc = ();

(my $otab = $iobo) =~ s/\.obo/.all.terms/g; $otab =~ s/^.*\///g;
(my $orel = $iobo) =~ s/\.obo/.relations/g; $orel =~ s/^.*\///g;
(my $ochi = $iobo) =~ s/\.obo/.children/g; $ochi =~ s/^.*\///g;
(my $osim = $iobo) =~ s/^.*\///g;
if($icut >= 0) { $osim =~ s/\.obo/.semsim.ge$icut.dat/g; }
else { my $ctag = abs($icut); $osim =~ s/\.obo/.semsim.le$ctag.dat/g; }
if($islim) { $orel =~ s/\.relations/.slim.relations/g; }

my $pot = my $por = my $poc = 0;
unless(-e $otab) { open TAB, ">$otab"; print TAB "#ID\tDesc\n"; $pot = 1; }
unless(-e $orel) { open REL, ">$orel"; $por = 1; }
unless(-e $ochi) { open CHI, ">$ochi"; $poc = 1; }

foreach my $term (@ts_terms) {
    if($pot) { print TAB "$term\t$id_2desc->{$term}\n"; }

    my @anc;
    if($islim) {
        @anc = sort {$term_level{$a} <=> $term_level{$b}} grep { exists $slim_ids{$_} } $obog->all_successors($term); }
    else {
        @anc = sort {$term_level{$a} <=> $term_level{$b}} $obog->all_successors($term); }

    map { $term_anc{$term}{$_}++; } @anc;
    $term_nanc{$term} = scalar @anc;

    if($por) {
        if((scalar @anc) == 0) { print "\n\tNo anc: $term\t$id_2desc->{$term}"; }
        else {
            print REL join "\t", ($term, $id_2desc->{$term}, @anc); print REL "\n"; } }

    if($poc) {
        my @chi = $obog->predecessors($term);
        if((scalar @chi) == 0) { print "\nNo chi: $term\t$id_2desc->{$term}"; }
        else {
            print CHI join "\t", ($term, $id_2desc->{$term}, @chi); print CHI "\n"; } } }

if($pot) { close TAB; }
if($por) { close REL; }
if($poc) { close CHI; }


if($isemsim) {
    my @lgs1 = ();
    if($ilgs1) {
        open LG, "$ilgs1" or die "Can't open $ilgs1!";
        while (<LG>) {
            if($_ =~ /^#/) { next; }
            chomp($_); my @p = split '\t', $_;
            push(@lgs1, $p[0]); }
        close LG; }
    else { @lgs1 = @ts_terms; }

    my @lgs2 = ();
    if($ilgs2) {
        open LG, "$ilgs2" or die "Can't open $ilgs2!";
        while (<LG>) {
            if($_ =~ /^#/) { next; }
            chomp($_); my @p = split '\t', $_;
            push(@lgs2, $p[0]); }
        close LG; }
    elsif($ilgs1) {
        @lgs2 = @lgs1; }
    else { @lgs2 = @ts_terms; }

    $time = runtime(); print "\n$time: Calc semantic similarity b/w terms ...";
    my $npairs = 0;
    open SIM, ">$osim";
    for(my $i=0; $i<=$#lgs1; $i++) {
        my $t1 = $lgs1[$i];
        unless(exists $term_nanc{$t1}) { next; }
        if($term_nanc{$t1} < 3) { next; }

        my $jstr;
        if($ilgs1) {
            if($ilgs2) { $jstr = 0; }
            else { $jstr = $i+1; } }
        else { $jstr = $i+1; }

        for(my $j=$jstr; $j<=$#lgs2; $j++) {
            my $t2 = $lgs2[$j]; if($t1 eq $t2) { next; }
            unless(exists $term_nanc{$t2}) { next; }
            if($term_nanc{$t2} < 3) { next; }

            $npairs++; unless($npairs % 1000000) {
                $time = runtime(); print "\n\t$time: $npairs ..."; }

            my @com = grep { exists $term_anc{$t2}{$_} } keys %{$term_anc{$t1}};
            my %uni = ();
            @uni{keys %{$term_anc{$t1}}} = values %{$term_anc{$t1}};
            @uni{keys %{$term_anc{$t2}}} = values %{$term_anc{$t2}};

            if((scalar @com) < 3) { next; }
            #my $num = 0; map { $num += $term_nanc{$_} } @com;
            my $num = scalar @com; if($num == 0) { next; }
            #my $den = 0; map { $den += $term_nanc{$_} } keys %uni;
            my $den = scalar keys %uni;

            my $sim = ($num/$den);
            if($icut >= 0) {
                if($sim < $icut) { next; } }
            else {
                if($sim > abs($icut)) { next; } }

            print SIM "$t1\t$t2\t", sprintf("%6f\n", $sim); } }
    close SIM; }


print "\n\tNo. terms: ", scalar @ts_terms;
unless($ogmt) {
    $time = runtime(); print "\nDONE $time\n\n";
    exit; }

# Mapping IDs
my %id_map = ();
if($imap) {
    open MAP, "$imap" or die "Can't open $imap!";
    while (<MAP>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $id_map{$p[0]} = $p[1]; }
    close MAP; }

# Parsing annotation file
$time = runtime(); print "\n$time: Parsing gene annotations ...";
my %term_genes = (); my %missing_terms = ();

open IANN, "$iann" or die "Can't open $iann!";
if($iann =~ /\.gmt$/) {
    while (<IANN>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;

        my $term = shift @p;
        unless((exists $id_2desc->{$term}) || (exists $alt2std_id->{$term})) {
            $missing_terms{$term}++;
            next; }

        $p[0] =~ s/ \([0-9]*\)$//g; $p[0] = lc($p[0]);
        $p[0] =~ s/ /_/g; $p[0] =~ s/__*/_/g;
        $id_2desc->{$term} = $p[0];

        shift @p; shift @p;
        foreach my $g (@p) {
            my $gene = $g;
            if($imap) {
                if(exists $id_map{$g}) { $gene = $id_map{$g}; }
                else { next; } }

            $term_genes{$term}{$gene}++; } } }
else {
    while (<IANN>) {
        if($_ =~ /^#/) { next; }
        chomp($_);

        my @p = split '\t', $_;
        my $term = $p[0]; my $gene = $p[1];
        if($imap) {
            if(exists $id_map{$gene}) { $gene = $id_map{$gene}; }
            else { next; } }

        if(exists $id_2desc->{$term}) {
            $term_genes{$term}{$gene}++; }
        elsif(exists $alt2std_id->{$term}) {
            $term_genes{$alt2std_id->{$term}}{$gene}++; }
        else { $missing_terms{$term}++; } } }
close IANN;

$time = runtime(); print "\n$time: Missing terms ...";
foreach (sort keys %missing_terms) {
    print "\n\t$_"; }

# Propagating gene annotations
$time = runtime(); print "\n$time: Propagating annotations ...";

my $num_terms = 0; my %final_genes = (); my @s;
my %prop_term_genes = ();
my %gene_slim = ();

open GMT, ">$ogmt";

foreach my $term (@ts_terms) {
    if($inoprop) {
        unless(exists $term_genes{$term}) { next; }
        print GMT "$term\t$id_2desc->{$term}"; $num_terms++;

        my $ngenes = scalar keys %{$term_genes{$term}};
        print GMT " ($ngenes)";

        foreach my $gene (keys %{$term_genes{$term}}) {
            $final_genes{$gene}++;
            print GMT "\t$gene"; }
        print GMT "\n"; }
    else {
        if(exists $term_genes{$term}) {
            foreach my $gene (keys %{$term_genes{$term}}) {
                $prop_term_genes{$term}{$gene}++; } }

        unless(exists $prop_term_genes{$term}) { next; }
        print GMT "$term\t$id_2desc->{$term}"; $num_terms++;

        # Assigning slim terms if this option is selected
        if($islim) {
            @s = grep { exists $slim_ids{$_} } $obog->all_successors($term);
            if(exists $slim_ids{$term}) { push(@s, $term); }
            if(scalar @s == 0) { print "\n\tNo slim: $term\t$id_2desc->{$term}"; }
            else {
                foreach my $t (@s) {
                    $t = $id_2desc->{$t}; }
                print GMT " ", join '|', @s; } }

        my $ngenes = scalar keys %{$prop_term_genes{$term}};
        print GMT " ($ngenes)";

        foreach my $gene (keys %{$prop_term_genes{$term}}) {
            $final_genes{$gene}++;
            print GMT "\t$gene";

            if($islim) {
                foreach my $t (@s) {
                    $gene_slim{$gene}{$t}++; } }

            my @parents = $obog->successors($term);
            foreach my $pt (@parents) {
                $prop_term_genes{$pt}{$gene}++; } }

        print GMT "\n"; } }
close GMT;

print "\n\n\tNo. Terms: $num_terms\n\tNo. Genes: ", scalar keys %final_genes;


if($islim) {
    my $ogs = $iann;
    if($iann =~ /\.(gmt|txt)$/) { $ogs =~ s/\.(gmt|txt)$/\.gene-slimcount.txt/g; }
    else { $ogs .= '.gene-slimcount.txt'; }

    open OGS, ">$ogs";
    foreach my $g (keys %gene_slim) {
        print OGS "$g\t", scalar keys %{$gene_slim{$g}}, "\n"; }
    close OGS; }


$time = runtime(); print "\nDONE $time\n\n";


# Subroutines
# ===========

sub parse_obo {
    my $file = shift;
    my $filt = shift;
    my $in = 0; my (@q, $term, $name, $rel, $parent);

    my $dag = Graph->new;
    my %id2desc = (); my %alt2id = (); # my %desc_2id = ();

    my $c = 0;

    open OBO, "$file" or die "Can't open $file!";
    LINE: while (<OBO>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @q = split ' ', $_;

        if ($_ =~ /^$/) { next LINE; }
        elsif ($q[0] eq '[Term]') { $in = 1; }
        elsif ($q[0] eq '[Typedef]') { $in = 0; }
        elsif ($in and ($q[0] eq 'id:')) {
            $term = $q[1];
            if($igsfilt) {
                unless(exists $filt->{$term}) {
                    $in = 0; } } }
        elsif ($in and ($q[0] eq 'name:')) {
            shift(@q); $name = lc(join ' ', @q);
            $name =~ s/[:;\/, ]/_/g; $name =~ s/[()'"]//g;
            $name =~ s/_-_/_/g; $name =~ s/__*/_/g;
            $id2desc{$term} = $name; }
        elsif ($in and ($q[0] eq 'alt_id:')) {
            # push(@{$alt2id{$q[1]}}, $term);
            $alt2id{$q[1]} = $term; }
        elsif ($in and ($q[0] eq 'is_a:')) {
            $parent = splice(@q, 0, 2);
            if($igsfilt) {
                unless(exists $filt->{$parent}) { next LINE; } }
            $dag->add_edge($term, $parent); }
        elsif ($in and ($q[0] eq 'relationship:')) {
            $rel = $q[1]; $parent = splice(@q, 0, 3);
            if($igsfilt) {
                unless(exists $filt->{$parent}) { next LINE; } }
            if($rel eq 'has_part') { next LINE; }
            elsif($rel eq 'part_of') {
                $dag->add_edge($term, $parent); }
            elsif(($rel eq 'develops_from') or ($rel eq 'related_to')) {
                $dag->add_edge($term, $parent); }
            elsif($iregp and ($rel =~ /regulates/)) {
                    $dag->add_edge($term, $parent); } }
        elsif ($in and ($q[0] eq 'is_obsolete:')) {
            $dag->delete_vertex($term);
            delete $id2desc{$term}; } }
    close OBO;

    return ($dag, \%id2desc, \%alt2id);
}


__END__

=head1

Propagates annotations using the true-path-rule.

=head1 USAGE

./propagate_annotations.pl [--iann DIRECT_ANNOTATIONS] [--iobo ONTOLOGY_OBO]
[--imap GENE_IDMAP] [--ogmt PROPAGATED_ANNOTATIONS] [--help]

=head1 DESCRIPTION

This script takes as input (i) a set of direct annotations of genes to terms in
an ontology, and (ii) an ontology file describing the annotation terms and
relationships, and propagates gene annotations satifying the ontology based on
'is-a' and 'part-of' relationships between terms.

=head1 ARGUMENTS

=over 12

=item C<--iann>

File containing direct annotations of genes to terms: <TermID> <GeneID>
tab-delimited; one/line.

=item C<--iobo>

Ontology file in OBO format.

=item C<--islim>

(Optional) Slim terms in the ontology: <TermID> <Desc> tab-delimited; one/line.
If this file is provided, then, (i) the slim term corresponding to a
given term is printed alongside the term description in the GMT file, and (ii)
an additional output file containing mapping of all terms to their slim terms is
created.

=item C<--imap>

Gene ID mapping file: <GeneID_old> <GeneID_new> tab-delimited; one/line.

=item C<--ogmt>

(Optional) If provided, this file is populated with gene annotations propagated
based on the ontology in GMT format. If not, all the term-term relationships
from the OBO file are printed. If --islim is also provided, then these
relationships are restricted to term to slim-term relationships.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Feb 25

=cut


