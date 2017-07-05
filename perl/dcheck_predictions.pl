#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

sub parse_gmt;
sub parse_obo;
sub read_fcol;

my($help, $imat, $igmt, $iobo, $iposslim, $inegslim, $itag, $igmt2, $otab);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'igmt=s' => \$igmt,
          'iobo=s' => \$iobo,
      'iposslim=s' => \$iposslim,
      'inegslim=s' => \$inegslim,
          'itag=s' => \$itag,
         'igmt2=s' => \$igmt2,
          'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $time, @gsref);

open PR, "$imat" or die "Can't open $imat!";
chomp(my @pmat = <PR>); close PR;

my @pred_func = split '\t', shift @pmat; shift @pred_func;
my %gene_pred = ();
foreach (@pmat) {
    @p = split '\t', $_;
    for(my $j=1; $j<=$#p; $j++) {
        if(($p[$j] eq 'NA') or ($p[$j] eq '') or ($p[$j] == 0)) { next; }
        $gene_pred{$pred_func[$j-1]}{$p[0]} = $p[$j];
    }
}

@gsref = parse_gmt($igmt);
my %gs_posg = %{$gsref[0]};
my %gs_desc = %{$gsref[1]};
my %gs_size = %{$gsref[2]};
my %all_genes = %{$gsref[3]};

my %all_slimg = ();
foreach my $gs (keys %gs_posg) {
    foreach my $g (keys %{%gs_posg{$gs}}) {
        $all_slimg{$g}++; } }

my %gs_posslim = ();
if($iposslim) {
    %gs_posslim = read_colf($iposslim); }
else {
    foreach my $gs (keys %gs_posg) {
        $gs_posslim{$gs}++; } }

if($inegslim) { my %gs_negslim = read_fcol($inegslim); }

my @a; my %gs_anc = ();
foreach my $gs (keys %gs_posslim) {
    @a = $obog->all_successors($gs);
    @a = grep { exists $gs_posslim{$_} } @a;
    if(scalar @a == 0) { next; }
    foreach (@a) { $gs_anc{$t}{$_}++; }
}

my %gs_negg = (); my $apar;
foreach my $gs (keys %gs_posslim) {
    foreach my $g (keys %all_slimg) {
        if(exists $gs_posg{$gs}{$g}) { next; }

        $apar = 1;
        ANC: foreach my $ags (keys %{$gs_anc{$gs}}) {
            if(exists $gs_posg{$ags}{$g}) { $apar = 0; last ANC; } }
        if($apar == 0) { next; }

        $gs_negg{$gs}{$g}++;
    }
}

foreach my $gs (keys %gene_pred) {
    # ...
}


# Subroutines
# Assigns genes, description & size to gs
sub parse_gmt {
    my $gmt = shift;

    my %genes = (); my %desc = (); my %size = ();
    my (@q, $gs); my %justgenes = ();

    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        @q = split '\t', $_;
        $gs = shift(@q);
        ($desc{$gs} = shift @q) =~ s/ \([0-9]*\)//g;

        foreach my $g (@q) {
            $genes{$gs}{$g}++;
            $justgenes{$g}++; }

        $size{$gs} = scalar keys %{$genes{$gs}};
    }
    close GMT;

    return (\%genes, \%desc, \%size, \%justgenes);
}

sub read_fcol {
    my $file = shift;
    my %col = ();
    open FH, "$file" or die "Can't open $file!";
    while (<FH>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $col{$p[0]}++;
    }
    close FH;

    return %col;
}

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
        elsif ($in and ($q[0] eq 'id:')) { ($term = $q[1]) =~ s/:/_/g; }
        elsif ($in and ($p[0] eq 'namespace:')) {
            if(($term =~ /^GO:/) and ($line !~ /biological_process/)) {
                $in = 0; next LINE; } }
        elsif ($in and ($q[0] eq 'is_a:')) {
            ($parent = splice(@q, 0, 2)) =~ s/:/_/g;
            $dag->add_edge($term, $parent);
        }
        elsif ($in and ($q[0] eq 'relationship:')) {
            $rel = $q[1]; ($parent = splice(@q, 0, 3)) =~ s/:/_/g;

            if($rel eq 'has_part') { next LINE; }
            elsif(($rel eq 'part_of') or ($rel =~ /regulates/)) {
                $dag->add_edge($term, $parent); }
            elsif(($rel eq 'develops_from') or ($rel eq 'related_to')) {
                $dag->add_edge($term, $parent); }
        }
    }

    return $dag;
}

sub dcheck {
    my $predl = shift; # Original_Class -> Prediction_Score
    my $predt = shift; # File to write out dcheker table
    my $FH = shift; # File handle to write out evaluation summary

    open PL, "$predl"; chomp(my @pred=<PL>); close PL;
    open PT, ">$predt";

    my $tp = my $fp = my $tn = my $fn = 0;
    my $P = my $N = my $pr = my $rc = my $srr = 0;
    my $auc = my $auprc = my $wauprc = my $p10r = my $p20r = my $p50r = 0;
    my $idx = 0; my %pred_score = (); my %true_label = ();
    my $in10 = my $in20 = my $in50 = 1;

    print PT "#Label\tScore\tTP\tFP\tTN\tFN\tPR\tRC\n";

    foreach (@pred) {
        @p = split '\t', $_;
        if($p[0] == 1) { $P++; }
        else { $N++; }
        $true_label{$idx} = $p[0];
        $pred_score{$idx} = $p[1];
        $idx++;
    }

    $idx = 0;
    foreach (sort {$pred_score{$b} <=> $pred_score{$a}} keys %pred_score) {
        if($true_label{$_} == 1) {
            $tp++; $rc = $tp / $P;

            if(($rc >= 0.10) and $in10) {
                $p10r = $tp / ($tp + $fp); $in10 = 0; }
            if(($rc >= 0.20) and $in20) {
                $p20r = $tp / ($tp + $fp); $in20 = 0; }
            if(($rc >= 0.50) and $in50) {
                $p50r = $tp / ($tp + $fp); $in50 = 0; }

            $srr += 1/($idx+1);
            $wauprc += $tp / ($tp + $fp) / ($idx + 1);
            $auprc += $tp / ($tp + $fp); }
        else {
            $fp++;
            $auc += $tp; }

        $idx++;
        $tn = $N - $fp; $fn = $P - $tp;
        $pr = $tp / ($tp + $fp); $rc = $tp / ($tp + $fn);

        print PT "$true_label{$_}\t$pred_score{$_}\t$tp\t$fp\t$tn\t$fn\t$pr\t$rc\n";
    }

    if(($tp == 0) or ($fp == 0)) {
        print "warning: Too few +ve true labels or -ve true labels\n";
        $auc = 0; $auprc = 0; $wauprc = 0;
    }
    else {
        $auc = $auc / $tp / $fp;
        $auprc = $auprc / $tp;
        $wauprc = $wauprc / $srr;
    }

    $predl =~ s/\.pred\.labels//g;
    print $FH "$predl\t$ppar\t$P\t$N\t", sprintf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
        $auc, $p10r, $p20r, $p50r, $auprc, $wauprc), "\n";

    close PL; close PT;

    `rm -f $pred_labels`;
}


__END__

=head1

<Brief_Desc>

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

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

