#!/usr/bin/perl
use strict;
use warnings;
use PDF;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);

sub gs2genes;
sub genes2fam;
sub get_fdr;

my $iqval = 0.01;
my($help, $time, $igmt1, $igmt2, $ifam, $ispc, $otab);
# igmt1: genesets in Species 1 .gmt
# igmt2: genesets in Species 2 .gmt
# ifam: gene-family information: <FamID> <s1g1> <s1g2> <s2g1> <s2g3> <s2g4>
# ispc: species pair, e.g. hsap:mmus
# otab: geneset pairs across species and their enrichment statistics

pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
         'igmt1=s' => \$igmt1,
         'igmt2=s' => \$igmt2,
          'ifam=s' => \$ifam,
          # 'ispc=s' => \$ispc,
          'iqval=s' => \$iqval,
          'otab=s' => \$otab ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

print "\nFile:\n\t*$igmt1*\n\t*$igmt2*\n\t*$ifam*\n";

open HH, ">$otab";


# Creating gene -> family hash
$time = runtime(); print "\n$time\tIndexing gene-family information ...";

# my ($s1, $s2) = split ":", $ispc;
# my %gene_fam = genes2fam($ifam, $s1, $s2);
my %gene_fam = genes2fam($ifam);

my $univ = scalar values %gene_fam;


# Creating gs -> fam hashes
$time = runtime(); print "\n$time\tIndexing geneset-gene information ...";

my ($gs1_genes, $gs1_desc) = genegs2($igmt1);
my ($gs1_gfam, $gs1_size) = gs_genes2fam($gs1_genes, \%gene_fam);

my ($gs2_genes, $gs2_desc) = genegs2($igmt2);
my ($gs2_gfam, $gs2_size) = gs_genes2fam($gs2_genes, \%gene_fam);


# Calculating overlap between all pairs of genesets between s1 and s2
$time = runtime(); print "\n$time\tCalculating geneset-overlaps ...";

my %gs_ncom = (); my %gs_pval = (); my %gs_lor = ();
my ($common, $tag, $lor); 
foreach my $gs1 (keys %gs1_gfam) {
    foreach my $gs2 (keys %gs2_gfam) {
        if($gs1 eq $gs2) { next; }

        $common = 0;
        foreach my $fam (keys %{$gs1_gfam{$gs1}}) {
            if(exists ${$gs2_gfam{$gs2}}{$fam}) {
                $common++;
            }
        }

        if($common < 3) { next; }
        $tag = $gs1.'__'.$gs2;

        $gs_ncom{$tag} = $common;

        $gs_pval{$tag} = hypergeometric_tail($univ, $gs1_size{$gs1},
            $gs2_size{$gs2}, $common);

        $lor = 0;   # Log-odds ratio
        unless($common == 0) {
            $lor = log(($common/$gs1_size{$gs1})/($gs2_size{$gs2}/$univ))/log(2); }
        $gs_lor{$tag} = $lor;
    }
}

# Performing multiple-hypothesis correction
$time = runtime(); print "\n$time\tPerforming multiple-hypothesis correction ...";

my %gs_fdr = get_fdr(\%gs_pval);

# Printing results
$time = runtime(); print "\n$time\tPrinting results ...";

my @p;
foreach my $sets (sort {$gs_fdr{$a} <=> $gs_fdr{$b}} keys %gs_fdr) {
    if($gs_fdr{$sets} > $iqval) { last; }

    my ($gs1, $gs2) = split "__", $sets;
    print HH "$gs1\t$gs1_desc{$gs1}\t$gs1_size{$gs1}\t";
    print HH "$gs2\t$gs2_desc{$gs2}\t$gs2_size{$gs2}\t";
    print HH "$gs_ncom{$sets}\t$gs_pval{$sets}\t";
    print HH "$gs_lor{$sets}\t$gs_fdr{$sets}\n"; }

close HH;

$time = runtime(); print "\n$time\tDONE\n\n";


# Assigns to family to genes
sub genes2fam {
    my $file = shift;

    open FH, "$file" or die "Can't open $file!";
        chomp(my @ortho=<FH>); close FH;

    my (@fl, %gene_fam, $fam, $o, $g);
    foreach (@ortho) {
        @fl = split '\t', $_;
        $fam = shift(@fl);
        foreach my $h (@fl) {
            ($o, $g) = split /\|/, $h;
            $gene_fam{$g} = $fam;
        }
    }

    return %gene_fam;
}

# Assigns genes and description to gs
sub gs2genes {
    my $file = shift;
    open FH, "$file" or die "Can't open $file!";
        chomp(my @gset=<FH>); close FH;

    my (@fl, %gs_genes, %gs_desc, $gs);
    foreach (@gset) {
        @fl = split '\t', $_;
        $gs = shift(@fl);
        # $fl[0] =~ s/ \([0-9]*\)//g;
        ($gs_desc{$gs} = shift(@fl)) =~ s/ \([0-9]*\)$//g;

        foreach my $g (@fl) {
            $gs_genes{$gs}{$g}++; } }

    return (\%gs_genes, \%gs_desc); }

# Maps genes to their gene-family IDs / meta-genes
sub gs_genes2fam {
    my $gs_genes = shift;
    my $gene_fam = shift;

    my %gs_gfam = ();
    my %gs_size = ();
    
    foreach my $gs (keys %{$gs_genes}) {
        foreach my $g (keys %{$gs_genes->{$gs}}) {
            unless(exists $gene_fam->{$g}) { next; }
            $gs_gfam{$gs}{$gene_fam->{$g}}++; }

        $gs_size{$gs} = scalar keys %{$gs_gfam{$gs}}; }

    return (\%gs_gfam, \%gs_size); }

# Performs multiple-hypothesis correction using BH method
sub get_fdr {
    my $pvalues_ref = shift;
    
    my $min_pval = 1; my $num_tests = 0;
    foreach my $sets (keys %{$pvalues_ref}) {
        $num_tests++;
        if(($pvalues_ref->{$sets} != 0) 
                and ($pvalues_ref->{$sets} < $min_pval)) {
            $min_pval = $pvalues_ref->{$sets};
        }
    }

    my %fdr = (); my $rank = $num_tests;
    foreach my $sets (sort {$pvalues_ref->{$b} <=> $pvalues_ref->{$a}}
        keys %{$pvalues_ref}) {
        if($pvalues_ref->{$sets} == 0) {
            $pvalues_ref->{$sets} = 0.1*sprintf("%.6g", $min_pval);
        }

        $fdr{$sets} = (($pvalues_ref->{$sets} * $num_tests) / $rank);
        if($fdr{$sets} > 1) { $fdr{$sets} = 1; }

        $rank--;
    }

    return %fdr;
}


__END__

=head1 NAME

code.pl
    - <Brief_Desc>

=head1 USAGE

./code.pl [--i INPUT_FILE] [--p OPTION] [--o OUTPUT_FILE] [--help]

=head1 DESCRIPTION

This script ...

=head1 ARGUMENTS

=over 12

=item C<--i>

Input file

=item C<--p>

Option

=item C<--o>

Output file

=item C<--help>

prints this documentation

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

=cut


