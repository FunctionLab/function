#!/usr/bin/perl

# ================================================
# Name : run_connectivity.pl
# Purpose : Run connectivity.pl in various settings
# Created : 03-02-2014
# Last Modified : Sat 22 Feb 2014 02:02:01 PM EST
# Author(s) : Arjun Krishnan
# ================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my ($help, $inet, $igmt, $islim, $igfilt, $igsfilt, $time, @p);
my $ineg = 0; my $icvk = 5; my $icvn = 10; my $ipred = 0;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
          'inet=s' => \$inet,
          'igmt=s' => \$igmt,
         'islim=s' => \$islim,
        'igfilt=s' => \$igfilt,
       'igsfilt=s' => \$igsfilt,
          'ineg=s' => \$ineg,
           'ipred' => \$ipred,
          'icvn=i' => \$icvn,
          'icvk=i' => \$icvk ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);




__END__

=head1

Run connectivity.pl in various settings.

=head1 USAGE

run_connectivity.pl [--inet NETWORK_DAB] [--igmt GENESETS_GMT] [--igfilt
GENE_FILTER] [--islim GENESETS_SLIM] [--ineg NEG_FLAG] [--icvk
CROSS-VALIDATION_FOLDS] [--icvn CROSS-VALIDATION_REPEATS] [--help]

=head1 DESCRIPTION

This script runs connectivity.pl, generating examples for each geneset in an
input collection, setting genes in the geneset as positives and other genes –
based on the islim, ineg and igfilt options – as negatives.

=head1 ARGUMENTS

=over 12

=item C<--inet>

Network file in DAB format.

=item C<--igmt>

Geneset collection in GMT format.

=item C<--igfilt>

(Optional) List of genes to include.

=item C<--islim>

(Optional) Geneset-slim associations that need to be used to set negatives.

=item C<--ineg>

(Optional) Flag to determine if negatives will be used in prediction. '0' if
–ves are to be used only for evaluation. '1' if they are to be used for
prediction too. Default '1'.

=item C<--iprior>

(Optional) Used to decide the number of negatives. Default 0.05.

=item C<--icvk>

(Optional) Cross-validation fold. Default 5.

=item C<--icvn>

(Optional) Cross-validation repeats. Default 10.

=item C<--ipred>

(Optional) Flag to only do prediction based on labelled examples without any CV.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2014 Feb 03

=cut

