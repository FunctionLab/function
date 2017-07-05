#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Graph;
use Time::SoFar qw(runtime);
#use Data::Dumper;

sub parse_obo;
sub gs2genes;
sub print_prop_ann;

my ($help, $igmt, $iobo, $irel, $igsincl, $igsfilt, $igfilt, $imaxg, $itinfo, $inonorm, $igenesim, $iexp, $ofile);
my $iming = 10; # my $in_glist = '';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<2);
GetOptions( 'help' => \$help,
          'iobo=s' => \$iobo,
          'irel=s' => \$irel,
       'igsincl=s' => \$igsincl,
       'igsfilt=s' => \$igsfilt,
        'igfilt=s' => \$igfilt,
          'igmt=s' => \$igmt,
        'itinfo=s' => \$itinfo,
         'inonorm' => \$inonorm,
      'igenesim=s' => \$igenesim,
            'iexp' => \$iexp,
         'iming=i' => \$iming,
         'imaxg=i' => \$imaxg  ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $time, @ref, $count);

# Parse OBO
$time = runtime(); print "\n$time: Parsing the OBO, GMT files ...";

my $obog;
if($igsincl) {
    $obog = parse_obo($iobo, $igsincl); }
else {
    $obog = parse_obo($iobo); }

# Parse GMT
my ($gene_terms, $term_desc, $term_size) = gs2genes($igmt);

# Term-info
$time = runtime(); print "\n$time: Calculating term-info ...";

## Calculate maximum term-info
my @term_array = sort {$term_size->{$a} <=> $term_size->{$b}} keys %$term_size;
my $tot_genes = scalar keys %$gene_terms;
my $maxi = (-1)*log(1/$tot_genes)/log(2);
#my $maxi = (-1)*log($iming/$tot_genes)/log(2);
#my $maxi = (-1)*log($term_size->{$term_array[0]}/$tot_genes)/log(2);
#my $mini = (-1)*log($term_size->{$term_array[$#term_array]}/$tot_genes)/log(2);
unless($imaxg) { $imaxg = $term_size->{$term_array[$#term_array]}; }


print "\n\tNo. genes: $tot_genes\n\tNo. terms: ", scalar @term_array;
print "\n\tMax term-info: $maxi";


my %term_succ = ();
foreach my $t (@term_array) {
    push(@{$term_succ{$t}}, $obog->all_successors($t)); }


## Reading/Printing term-info
$ofile = $igmt; $ofile =~ s/^.*\///g;
$ofile =~ s/\.gmt/\.term-info\.txt/g;
#print "\n\t$iming term-info: $ofile";

my %term_info = ();
if($itinfo) {
    open TINF, "$itinfo" or die "Can't open $itinfo!";
    while (<TINF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $term_info{$p[0]} = $p[1]; }
    close TINF;
} elsif($irel) {
    my %counts = my %desc = ();
    open REL, "$irel" or die "Can't open $irel!";
    while (<REL>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $desc{$p[0]}++;
        foreach my $anc (@p[2..$#p]) {
            $counts{$anc}++; } }
    close REL;

    foreach my $d (keys %desc) {
        if(exists $counts{$d}) { next; }
        $counts{$d}++; }

    my @sort_counts = sort {$b <=> $a} values %counts;
    my $max_counts = $sort_counts[0];

    foreach my $term (keys %counts) {
        $term_info{$term} = (-1)*log($counts{$term}/$max_counts)/log(2); }

    open TINF, ">$ofile";
    print TINF "#ID\tInfo\tSize\tDesc\n";
    foreach my $t (sort {$term_info{$b} <=> $term_info{$a}} keys %term_info) {
        print TINF "$t\t", sprintf("%.6f", $term_info{$t}), "\n"; }
    #print TINF "$t\t", sprintf("%.6f", $term_info{$t}), "\t$term_size->{$t}\t$term_desc->{$t}\n"; }
    close TINF;
} elsif(-e $ofile) {
    print "\n\t$ofile already exists!";
    open TINF, "$ofile" or die "Can't open $ofile!";
    while (<TINF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $term_info{$p[0]} = $p[1]; }
    close TINF;
} else {
    # Generate term-array & calculate normalized term-info
    foreach my $t (@term_array) {
        $term_info{$t} = (-1)*log($term_size->{$t}/$tot_genes)/log(2);
        unless($inonorm) { $term_info{$t} /= $maxi; }; }

    #%term_info = equalizetermswidentanc(\%term_info);

    open TINF, ">$ofile";
    print TINF "#ID\tInfo\tSize\tDesc\n";
    foreach my $t (sort {$term_info{$b} <=> $term_info{$a}} keys %term_info) {
        print TINF "$t\t", sprintf("%.6f", $term_info{$t}), "\t$term_size->{$t}\t$term_desc->{$t}\n"; }
    close TINF;
}


my %gs_filt = ();
if($igsfilt) {
    open GSFILT, "$igsfilt" or die "Can't open $igsfilt!";
    while (<GSFILT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $gs_filt{$p[0]}++; }
    close GSFILT; }
else {
    foreach my $gs (keys %$term_size) {
        if(($term_size->{$gs} < $iming) or ($term_size->{$gs} > $imaxg)) { next; }
        $gs_filt{$gs}++; } }


# Term-similarity
$time = runtime(); print "\n$time: Calculating term-term sim ...";
$ofile =~ s/term-info\.txt/term-sim\.dat/g;
#print "\n\tterm-sim: $ofile"; exit;

my %term_sim = (); my ($t1, $t2, @a, %anc1, %anc2);

if(-e $ofile) {
    print "\n\t$ofile already exists!";

    if($igenesim) {
        open TSIM, "$ofile" or die "Can't open $ofile!";
        while (<TSIM>) {
            if($_ =~ /^#/) { next; }
            chomp($_); @p = split '\t', $_;
            $term_sim{join '__', sort ($p[0], $p[1])} = $p[2];
        }
        close TSIM;
    }
}
else {
    @term_array = grep { exists $gs_filt{$_}; } @term_array;
    open TSIM, ">$ofile";

    # Calculate term-term similarity
    $count = 0;
    for (my $i=0; $i<$#term_array; $i++) {
        $t1 = $term_array[$i];
        $term_sim{join '__', ($t1, $t1)} = $term_info{$t1};

        #@a = $obog->all_successors($t1);
        @a = grep { exists $term_info{$_} } @{$term_succ{$t1}}; push(@a, $t1);
        %anc1 = (); foreach (@a) { $anc1{$_}++; }

        $count++;
        for (my $j=($i+1); $j<=$#term_array; $j++) {
            $t2 = $term_array[$j];

            #@a = $obog->all_successors($t2);
            @a = grep { exists $term_info{$_} } @{$term_succ{$t2}}; push(@a, $t2);
            %anc2 = (); foreach (@a) { $anc2{$_}++; }

            my @com = grep { exists $anc1{$_} } keys %anc2;
            if((scalar @com) == 0) { next; }

            @a = sort {$term_info{$b} <=> $term_info{$a}} @com;
            $maxi = $term_info{$a[0]};

            my %uni = ();
            @uni{keys %anc1} = values %anc1;
            @uni{keys %anc2} = values %anc2;

            foreach my $t (keys %uni) {
                unless(exists $term_info{$t}) {
                    print "\nno term info for $t\n"; exit; } }
            my $num = 0; map { $num += $term_info{$_} } @com;
            my $den = 0; map { $den += $term_info{$_} } keys %uni;

            if($maxi > 0) {
                print TSIM "$t1\t$t2\t", sprintf("%.6f", $maxi),"\n"; }
        }

        if(($count % 100)==0) {
            $time = runtime();
            print "\n  $time: $count\t\t$t1 $t2 $maxi";
        }
    }

    $t1 = $term_array[$#term_array];
    $term_sim{join '__', ($t1, $t1)} = $term_info{$t1};

    close TSIM;
}

unless($igenesim) {
    $time = runtime(); print "\nDONE $time\n\n"; exit; }

if($iexp) {
    $igmt =~ s/\.closed\.gmt/\.EXP\.closed\.gmt/g;
    print "\n\t$igmt ... ";
    if(-e $igmt) { ($gene_terms, $term_desc, $term_size) = gs2genes($igmt); }
    else { print "does not exist!"; }
}


my %gene_filt = ();
if($igfilt) {
    open GFIL, "$igfilt" or die "Can't open $igfilt!";
    while (<GFIL>) {
        if($_ =~ /^#/) { next; }
        chomp($_); $gene_filt{$_}++; }
    close GFIL; }
else {
    foreach my $g (keys %$gene_terms) {
        my @at = grep { exists $gs_filt{$_}; } keys %{$gene_terms->{$g}};
        if((scalar @at) == 0) {
            delete $gene_terms->{$g}; next; }

        my @notat = grep { not exists $gs_filt{$_}; } keys %{$gene_terms->{$g}};
        if((scalar @notat) == 0) { next; }
        map { delete $gene_terms->{$g}->{$_}; } @notat; } }

my @gene_array = keys %$gene_terms;
if($igfilt) { @gene_array = grep { exists $gene_filt{$_}; } @gene_array; }
print "\n\tNo. genes: ", scalar @gene_array;


# Topological sorting and propagation
$time = runtime(); print "\n$time: Topologically sorting all the terms ...";
my @ts_terms = $obog->topological_sort;
my %term_order = (); $count = 0;
foreach (@ts_terms) { $term_order{$_} = $count; }

# Retain only leaf-terms for each gene
$time = runtime(); print "\n$time: Retaining only leaf-terms for each gene ...";

# my @des; $count = 0;
my @suc; $count = 0;
foreach my $g (@gene_array) {
    TERM: foreach my $t (sort {$term_order{$a} <=> $term_order{$b}} keys %{$gene_terms->{$g}}) {
        @suc = grep { exists $gene_terms->{$g}->{$_} } @{$term_succ{$t}};
        foreach(@suc) { delete $gene_terms->{$g}->{$_}; } }

    $count++; unless($count % 1000) {
        $time = runtime(); print "\n\t$time: $count ..."; } }


# Calculate gene-gene similarity
$time = runtime(); print "\n$time: Calculating gene-gene sim ...\n";

($ofile = $igmt) =~ s/^.*\///g;
$ofile =~ s/\.gmt/\.$iming-$imaxg\.$igenesim\.gene-sim\.dab/g;
if($igsfilt) { $igsfilt =~ s/^\.*\///g; $ofile =~ s/\.gene-sim/\.$igsfilt\.gene-sim/g; }
#print "\n\t*$ofile*"; exit;

my ($g1, $g2, $genesim, $ttag, $gdat, $gdab); $count = 0;
my (@comt, %unit, $num, $den); my $nnon0 = 0;
open GSIM, ">gene.dat";
for(my $i=0; $i<$#gene_array; $i++) {
    $g1 = $gene_array[$i];
    $count++;

    for(my $j=($i+1); $j<=$#gene_array; $j++) {
        $g2 = $gene_array[$j];
        $genesim = 0;

        if($igenesim eq 'max') {
            foreach my $t1 (keys %{$gene_terms->{$g1}}) {
                foreach my $t2 (keys %{$gene_terms->{$g2}}) {
                    $ttag = join '__', sort ($t1, $t2);
                    unless(exists $term_sim{$ttag}) { next; }
                    if($term_sim{$ttag} > $genesim) {
                        $genesim = $term_sim{$ttag}; } } } }

        elsif($igenesim eq 'jac') {
            @comt = grep { exists $gene_terms->{$g2}->{$_} }
                keys %{$gene_terms->{$g1}};

            %unit = ();
            @unit{keys %{$gene_terms->{$g1}}} = values %{$gene_terms->{$g1}};
            @unit{keys %{$gene_terms->{$g2}}} = values %{$gene_terms->{$g2}};

            $num = 0; $den = 0;
            foreach (@comt) { $num += $term_info{$_}; }
            foreach (keys %unit) { $den += $term_info{$_}; }

            $genesim = ($num/$den); }

        if($genesim > 0) {
            $nnon0++;
            print GSIM "$g1\t$g2\t", sprintf("%.6f", $genesim), "\n"; }
    }
    
#     if($nnon0 > 2500000) {
#         $time = runtime(); print "\t$time: $count\t$nnon0\n";
#         close GSIM; exit;
#         `Dat2Dab -i gene.dat -o gene.dab`;
# 
#         if(-e $ofile) {
#             `cp $ofile tmp.dab`;
#             `Combiner -t dat -o $ofile tmp.dab gene.dab`;
#             `rm -f gene.dat gene.dab tmp.dab`; }
#         else { `mv gene.dab $ofile; rm -f gene.dat`; }
# 
#         $nnon0 = 0;
#         open GSIM, ">gene.dat"; }
}

close GSIM;
# if((-e 'gene.dat') and ($nnon0 > 0)) {
#     close GSIM;
#     `Dat2Dab -i gene.dat -o gene.dab`;
#     if(-e $ofile) {
#         `cp $ofile tmp.dab`;
#         `Combiner -t dat -o $ofile tmp.dab gene.dab`;
#         `rm -f gene.dat gene.dab tmp.dab`; }
#     else { `mv gene.dab $ofile; rm -f gene.dat`; } }
    #`cp $ofile tmp.dab`;
    #`Combiner -t dat -o $ofile tmp.dab gene.dab`;
    #`rm -f gene.dat gene.dab tmp.dab`; }

$time = runtime(); print "\nDONE $time\n\n";

# Subroutines
# ===========
# Parse OBO to return graph object and a hash of all terms
sub parse_obo {
    my $obofile = shift;
    my $gslist = shift;
    open OBO, "$obofile" or die "Can't open $obofile!";
        chomp(my @obo=<OBO>); close OBO;
        
    my @p; my %sel_gs = ();
    if($gslist) {
        open GSL, "$gslist" or die "Can't open $gslist!";
        while(<GSL>) { chomp($_); @p = split '\t', $_; $sel_gs{$p[0]}++; }
        close GSL; }

    my $dag = Graph->new;

    my $inside = 0;
    my ($term, $name, $syn, $rel, $parent);
    LINE: foreach my $line (@obo) {
        if($line =~ /^#/) { next LINE; }
        @p = split(' ', $line);

        if ($line =~ /^$/) { next LINE; }
        elsif ($p[0] eq '[Term]') { $inside = 1; }
        elsif ($p[0] eq '[Typedef]') { $inside = 0; }
        elsif ($inside and ($p[0] eq 'id:')) {
            if($gslist) {
                if(exists $sel_gs{$p[1]}) { $term = $p[1]; }
                else { $inside = 0; } }
            else { $term = $p[1]; }
        }
        elsif ($inside and ($p[0] eq 'is_a:')) {
            $parent = splice(@p, 0, 2);
            if($gslist) { unless(exists $sel_gs{$parent}) { next LINE; } }
            $dag->add_edge($term, $parent);
        }
        elsif ($inside and ($p[0] eq 'relationship:')) {
            $rel = $p[1]; $parent = splice(@p, 0, 3);
            if($gslist) { unless(exists $sel_gs{$parent}) { next LINE; } }

            if($rel eq 'has_part') { next LINE; }
            elsif(($rel eq 'part_of') or ($rel =~ /regulates/)) {
                $dag->add_edge($term, $parent); }
            elsif(($rel eq 'develops_from') or ($rel eq 'related_to')) {
                $dag->add_edge($term, $parent); }
        }
        elsif ($inside and ($p[0] eq 'is_obsolete:')) {
            $dag->delete_vertex($term); $inside = 0; }
    }

    return $dag;
}

# Get terms in the ontology that have identical ancestors
sub equalizetermswidentanc {
    my $tiref = shift;

    my @tary = keys %$tiref;
    my $ttg = Graph::Undirected->new;

    # Link terms with identical ancestors
    for (my $i=0; $i<$#tary; $i++) {
        my $t1 = $tary[$i];
        unless($i % 1000) { print "\n\t\t$i ... $t1"; }

        my %anc1 = (); foreach (@{$term_succ{$t1}}) { $anc1{$_}++; }
        my $na1 = scalar @{$term_succ{$t1}};

        for (my $j=($i+1); $j<=$#tary; $j++) {
            my $t2 = $tary[$j];
            my $na2 = scalar @{$term_succ{$t2}};
            if($na2 != $na1) { next; }

            my %anc2 = (); foreach (@{$term_succ{$t2}}) { $anc2{$_}++; }
            my @com = grep { exists $anc1{$_} } keys %anc2;
            if((scalar @com) < (scalar @{$term_succ{$t2}})) { next; }

            $ttg->add_edge($t1, $t2); } }

    my @tcc = $ttg->connected_components();

    my %term_info = ();
    foreach my $c (@tcc) {
        my @ainfo = (); map { push(@ainfo, $tiref->{$_}); } @$c;
        @ainfo = sort {$a <=> $b} @ainfo;
        my $mininfo = $ainfo[0];

        foreach my $term (@$c) {
            $term_info{$term} = $mininfo; } }

    return %term_info; }

# Assigns genes and description to gs
sub gs2genes {
    my $gmt = shift;
    my (@p, %genes, %desc, %size, $gs);

    open GMT, "$gmt" or die "Can't open $gmt!";
    GSET: while (<GMT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $p[1] =~ s/ \([0-9]*\)$//g;

        $gs = shift(@p);
        $desc{$gs} = shift(@p);

        #if(scalar(@p) < $iming) { next GSET; }
        #elsif(scalar(@p) > $imaxg) { next GSET; }
        foreach my $g (@p) {
            $size{$gs}++; $genes{$g}{$gs}++; } }
    close GMT;

    return (\%genes, \%desc, \%size);
}

# Print propagated gene annotations
sub print_prop_ann {
    my $term_ref = shift;
    my $term_gene_ref = shift;
    my $dag = shift;
    my $fh = shift;

    my $num_terms = 0; my $num_genes = 0;
    my %final_genes = (); my %prop_term_genes = (); my  @parents;

    foreach my $term (@{$term_ref}) {
        if(exists $term_gene_ref->{$term}) {
            foreach my $gene (keys %{$term_gene_ref->{$term}}) {
                $prop_term_genes{$term}{$gene}++;
            }
        }

        unless(exists $prop_term_genes{$term}) { next; }
        $num_terms++;

        $num_genes = scalar(keys %{$prop_term_genes{$term}});

        foreach my $gene (keys %{$prop_term_genes{$term}}) {
            $final_genes{$gene}++;

            @parents = $dag->successors($term);
            foreach my $pt (@parents) {
                $prop_term_genes{$pt}{$gene}++;
            }
        }
    }

    $num_genes = scalar(keys %final_genes);
    print "\nNo. Terms: $num_terms\nNo. Genes: $num_genes\n\n";
}

__END__

=head1

Calculates semantic-similarity measures between GO term and gene pairs.

=head1 USAGE

./calculate_semantic-similarity.pl [--iobo ONTOLOGY_OBO] [--igsincl
SELECT_TERM-LIST] [--igmt DIRECT_ANNOTATIONS] [--help]

=head1 DESCRIPTION

This script takes (i) annotations of genes to terms in an ontology, &
(ii) ontology of terms and their relationships, to give semantic-similarity
measures (Resnik) between pairs of GO terms and genes.

=head1 ARGUMENTS

=over 12

=item C<--iobo>

Ontology file in OBO format.

=item C<--igsincl>

Select terms from the ontology, for e.g., those terms that belong to a
particular namespace (like 'biological_process' for GO). Format <TermID>
<TermDesc>

=item C<--igmt>

Gene annotations file in GMT format. Annotations should be propagated.
Currently, the code assumes that what is provided here is the set of annotations
based on all evidences (used to calcualte term information). For calculating
similarity between genes (i.e. with option --igenesim), if the iexp flag is
provided, the code can also look for another GMT file in the same folder that
contains annotations only based on experiments. If present, this second GMT is
used to restrict gene-term associations to only those that are reliable.

=item C<--igenesim>

(Optional) If given, calculate semantic similarity between genes. The options
are 'max' and 'jac'. If 'max' is chosen, then gene-sim between two genes
i and j = term-sim(argmax_(ti, tj) (term-sim(ti, tj))), where ti and tj are
terms annotating genes i and j. If 'jac' is chosen, then gene-sim =
weighted.Jaccard (Ti, Tj), where Ti and Tj are the sets of terms annotating
genes i and j.

=item C<--iexp>

(Optional) If given, then direct annotations only based on experimental
evidences is used for calculating similarity between genes.

=item C<--iming>

Size of the smallest term to be considered. (Optional; Default 10)

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Apr 20

=cut

