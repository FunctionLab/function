#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);
use Graph;

sub parse_gmt;
sub parse_obo;

my ($help, $igmt, $iobo, $oterm, $ogene);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'igmt=s' => \$igmt,
          'iobo=s' => \$iobo,
         'oterm=s' => \$oterm,
         'ogene=s' => \$ogene) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $time);

# Parsing GMT & OBO
$time = runtime(); print "\n$time: Parsing GMT & OBO files ...";

my ($gene_terms, $term_desc, $term_size) = parse_gmt($igmt);
my $obog = parse_obo($iobo);

# Calculate maximum term-info
$time = runtime(); print "\n$time: Calculating max-info ...";

my @term_array = sort {$term_size->{$a} <=> $term_size->{$b}} keys %$term_size;
my $tot_genes = scalar keys %{$gene_terms};
my $maxi = (-1)*log($term_size->{$term_array[0]}/$tot_genes);

# Generate term-array & calculate normalized term-info
$time = runtime(); print "\n$time: Calculating normlzd term-info ...";

open TINF, ">$oterm"; print TINF "#ID\tInfo\tSize\tDesc\n";
my %term_info = ();

foreach my $t (@term_array) {
    $term_info{$t} = (-1)*log($term_size->{$t}/$tot_genes)/$maxi;
    
    print TINF "$t\t", sprintf("%.3f", $term_info{$t});
    print TINF "\t$term_size->{$t}\t$term_desc->{$t}\n";
}
close TINF;

# Calculate term-term similarity
$time = runtime(); print "\n$time: Calculating term-term semsim ...";

my %terms_semsim = (); my ($t1, $t2, @a, %anc1); my $count = 0;
for (my $i=0; $i<$#term_array; $i++) {
    $t1 = $term_array[$i];
    $terms_semsim{join '__', ($t1, $t1)} = $term_info{$t1};

    @a = $obog->all_successors($t1);
    @a = grep { exists $term_info{$_} } @a;
    %anc1 = (); $anc1{$t1}++; foreach (@a) { $anc1{$_}++; }

    $count++;
    for (my $j=($i+1); $j<=$#term_array; $j++) {
        $t2 = $term_array[$j];

        @a = $obog->all_successors($t2); push(@a, $t2);
        @a = grep { exists $anc1{$_} } @a;
        @a = sort {$term_info{$b} <=> $term_info{$a}} @a;

        $terms_semsim{join '__', sort ($t1, $t2)} = $term_info{$a[0]};
    }

    unless($count % 100) {
        $time = runtime(); print "\n\t$time: $count";
    }
}

$t1 = $term_array[$#term_array];
$terms_semsim{join '__', ($t1, $t1)} = $term_info{$t1};

# Calculate gene multifunctionality
$time = runtime(); print "\n$time: Calculating gene multifunctionality scores ...";

open MSIM, ">$ogene";
print MSIM "#Gene\tNo.Leaf.Terms\tGP.Score\tSum.Info\tSum.SemDis\tMean.SemDis\tSemDis.Score\n";

my ($nterms, $wght_nterms, $semdis, $sum_info, $sum_semdis, $mean_semdis, $sd_semdis, $semdis_score);
my (@des, @leaves, $npairs, $ttag, $prev_mean, %terms2del); $count = 0;
my @all_semdis = (); my %remgene_nterms = (); my %remgene_npairs = ();
my %remgene_suminfo = (); my %remgene_sumsemdis = (); my %remgene_semdis = (); my %remgene_wnterms = ();
GENE: foreach my $g (keys %{$gene_terms}) {
    $count++; %terms2del = ();
    TERM: foreach my $t (keys %{$gene_terms->{$g}}) {
        @des = $obog->all_predecessors($t);
        DESC: foreach my $dt (@des) {
            if(exists $gene_terms->{$g}->{$dt}) {
                $terms2del{$t}++; next TERM; }
        }
    }

    @leaves = grep { not exists $terms2del{$_} } keys %{$gene_terms->{$g}};
    $nterms = scalar @leaves;
    
    $wght_nterms = 0; $npairs = 0; $sum_info = 0;
    $sum_semdis = 0; $mean_semdis = 0; $prev_mean = 0; $sd_semdis = 0;

    if($nterms == 1) {
        if($term_size->{$leaves[0]} == $tot_genes) { next GENE; }
        $wght_nterms = 1/($term_size->{$leaves[0]}*($tot_genes - $term_size->{$leaves[0]}));
        $wght_nterms = sprintf("%.5f", log($wght_nterms)/log(10));
        $sum_info = sprintf("%.5f", $term_info{$leaves[0]});
        $semdis_score = 0;
        print MSIM "$g\t$nterms\t$wght_nterms\t$sum_info\t$sum_semdis\t$mean_semdis\t$semdis_score\n";
        next GENE;
    }

    foreach my $t (@leaves) {
        $wght_nterms += 1/($term_size->{$t}*($tot_genes - $term_size->{$t}));
        $sum_info += $term_info{$t}; }
    $wght_nterms = sprintf("%.5f", log($wght_nterms)/log(10));
    $sum_info = sprintf("%.5f", $sum_info);

    for(my $i=0; $i<$#leaves; $i++) {
        for(my $j=($i+1); $j<=$#leaves; $j++) {
            $npairs++;
            $ttag = join '__', sort ($leaves[$i], $leaves[$j]);
            $semdis = 1 - $terms_semsim{$ttag}; push(@all_semdis, $semdis);
            $sum_semdis += $semdis;
            $mean_semdis = $prev_mean + ($semdis - $prev_mean)/$npairs;
            $sd_semdis += ($semdis - $mean_semdis)*($semdis - $prev_mean);
            $prev_mean = $mean_semdis;
        }
    }

    if($npairs == 1) {
        $remgene_nterms{$g} = $nterms;
        $remgene_npairs{$g} = $npairs;
        $remgene_semdis{$g} = $mean_semdis;
        $remgene_wnterms{$g} = $wght_nterms;
        $remgene_suminfo{$g} = $sum_info;
        $remgene_sumsemdis{$g} = $sum_semdis;
        next GENE; }
    
    if($sd_semdis == 0) { $semdis_score = 'NA'; }
    else {
        $sd_semdis = sqrt($sd_semdis/($npairs-1));
        $semdis_score = sprintf("%.5g", $mean_semdis/$sd_semdis);
    }
    
    $sum_semdis = sprintf("%.5g", $sum_semdis);
    $mean_semdis = sprintf("%.5g", $mean_semdis);
    $sd_semdis = sprintf("%.5g", $sd_semdis);

    print MSIM "$g\t$nterms\t$wght_nterms\t$sum_info\t$sum_semdis\t$mean_semdis\t$semdis_score\n";

    unless($count % 500) {
        $time = runtime(); print "\n\t$time: $count";
    }
}

$mean_semdis = 0; $prev_mean = 0; $sd_semdis = 0;
foreach (@all_semdis) {
    $mean_semdis = $prev_mean + ($_ - $prev_mean)/$npairs;
    $sd_semdis += ($_ - $mean_semdis)*($_ - $prev_mean);
    $prev_mean = $mean_semdis; }
$sd_semdis = sqrt($sd_semdis/(scalar @all_semdis - 1));

foreach my $g (keys %remgene_semdis) {
    $semdis_score = sprintf("%.5g",
        $remgene_semdis{$g}/$sd_semdis);
    # $wght_nterms = sprintf("%.5g", $remgene_wnterms{$g});
    $wght_nterms = $remgene_wnterms{$g};
    $sum_info = $remgene_suminfo{$g};
    $sum_semdis = sprintf("%.5g", $remgene_sumsemdis{$g});
    $semdis = sprintf("%.5g", $remgene_semdis{$g});

    print MSIM "$g\t$remgene_nterms{$g}\t$wght_nterms\t$sum_info\t";
    print MSIM "$sum_semdis\t$semdis\t$semdis_score\n";
}

close MSIM;

$time = runtime(); print "\n$time: DONE\n\n";

# Subroutines
# Parse OBO to return graph object and a hash of all terms
sub parse_obo {
    my $obofile = shift;
    open OBO, "$obofile" or die "Can't open $obofile!";
        chomp(my @obo=<OBO>); close OBO;
        
    my $dag = Graph->new;

    my $in = 0;
    my ($term, $name, $syn, $rel, $parent, @q);
    LINE: foreach my $line (@obo) {
        if($line =~ /^#/) { next LINE; }
        @q = split(' ', $line);

        if ($line =~ /^$/) { $in = 0; next LINE; }
        elsif ($q[0] eq '[Term]') { $in = 1; }
        elsif ($q[0] eq '[Typedef]') { $in = 0; }
        elsif ($in and ($q[0] eq 'id:')) { $term = $q[1]; }
        elsif ($in and ($q[0] eq 'namespace:')) {
            unless($line =~ /biological_process/) { $in = 0; } }
        elsif ($in and ($q[0] eq 'is_a:')) {
            $parent = splice(@q, 0, 2);
            $dag->add_edge($term, $parent);
        }
        elsif ($in and ($q[0] eq 'relationship:')) {
            $rel = $q[1]; $parent = splice(@q, 0, 3);

            if($rel eq 'has_part') { next LINE; }
            elsif(($rel eq 'part_of') or ($rel =~ /regulates/)) {
                $dag->add_edge($term, $parent); }
            elsif(($rel eq 'develops_from') or ($rel eq 'related_to')) {
                $dag->add_edge($term, $parent); }
        }
    }

    return $dag;
}

# Assigns genes, description & size to gs
sub parse_gmt {
    my $gmt = shift;

    my %genes = (); my %desc = (); my %size = ();
    my (@q, $gs); my $msize = 10;

    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        chomp($_); @q = split '\t', $_;
        $gs = shift(@q);
        $q[0] =~ s/ \([0-9]*\)//g;
        $desc{$gs} = shift(@q);

        if(scalar(@q) < $msize) { delete $desc{$gs}; next; }
        foreach my $g (@q) {
            $genes{$g}{$gs}++; }

        $size{$gs} = scalar @q;
    }
    close GMT;

    return (\%genes, \%desc, \%size);
}


__END__

=head1

Calculates a multifunctionality score for each of the annotated genes.

=head1 USAGE

calculate_multifunctionality.pl [--igmt ANNOTATIONS_GMT] [--iobo ONTOLOGY_OBO] [--otab
OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes in a complete annotations file in GMT format + an ontology file
in OBO format to calculate a multifunctionality score for each gene based on (i)
no. of associated GO terms, (ii) weighted count of associated GO terms, & (iii)
semantic similarity of associated GO terms.

=head1 ARGUMENTS

=over 12

=item C<--igmt>

Annotations in GMT format.

=item C<--iobo>

Ontology in OBO format

=item C<--oterm>

Output table of terms & their info scores.

=item C<--ogene>

Output table of genes & different measures.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 Jul 10

=cut

