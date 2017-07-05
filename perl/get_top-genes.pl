#!/usr/bin/perl

# ================================================
# Name : get_top-genes.pl
# Purpose : Given a MAT, get top-n genes per column and output a GMT
# Created : 07-02-2014
# Last Modified : Sun 09 Feb 2014 06:17:17 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, $imat, $ibot, $ogmt); my $inum = '10,25,50,100,250,500,1000';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'imat=s' => \$imat,
          'inum=s' => \$inum,
            'ibot' => \$ibot,
          'ogmt=s' => \$ogmt    ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my @num = split ',', $inum;


my %all_genes = my @col = my %gene_score = (); my $nrow = 0;
open MAT, "$imat" or die "Can't open $imat!";
while (<MAT>) {
    if($_ =~ /^#/) { next; }
    $nrow++;
    chomp($_); my @p = split '\t', $_;

    if($nrow == 1) { @col = @p; next; }
    $all_genes{$p[0]}++;
    for(my $j=1; $j<=$#p; $j++) {
        if($p[$j] =~ /^(nan|NaN|NA)$/) { next; }
        $gene_score{$col[$j]}{$p[0]} = $p[$j]; } }
close MAT;


open GMT, ">$ogmt";
print GMT "ALLG\tall_genes.", scalar keys %all_genes;
foreach my $g (keys %all_genes) {
    print GMT "\t$g"; }
print GMT "\n";

my $ncol = 0;
foreach my $c (sort keys %gene_score) {
    $ncol++;
    
    my @sorta = sort {$gene_score{$c}{$b} <=> $gene_score{$c}{$a}} keys %{$gene_score{$c}};
    foreach my $n (@num) {
        if($n > 0.75*($nrow-1)) { last; }

        my $id = 'COL'.sprintf("%02d", $ncol).'.TOP'.sprintf("%04d", $n);
        my $ds = $c.'.top'.$n;
        print GMT "$id\t$ds\t", join "\t", @sorta[0..($n-1)], "\n";

        unless($ibot) { next; }
        $id = 'COL'.sprintf("%02d", $ncol).'.BOT'.sprintf("%04d", $n);
        $ds = $c.'.bot'.$n;
        print GMT "$id\t$ds\t", join "\t", @sorta[($#sorta-$n+1)..$#sorta], "\n"; } }

close GMT;



__END__

=head1

Get top-n genes for each column in a matrix.

=head1 USAGE

get_top-genes.pl [--imat INPUT_MAT] [--ibot INCLUDE_BOTTOM] [--ogmt OUTPUT_GMT] [--help]

=head1 DESCRIPTION

This script takes in a matrix and output a GMT containing genesets corresponding
to the top-n gene per column for n=10,25,50,100,250,500,1000. Bottom genes will
also be included if option --ibot is provided.

=head1 ARGUMENTS

=over 12

=item C<--imat>

Input matrix with genes along the rows, with the first row containing column
names.

=item C<--ibot>

(Optional) Option to select bottom-n genes along with top genes.

=item C<--ogmt>

Output GMT file containing the top (& bottom) genesets per column.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2014 Feb 07

=cut

