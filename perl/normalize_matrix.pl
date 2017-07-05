#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::Distributions;

pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions('help' => \$help, 'i=s' => \$infile, 'p=s{,}' => \$option, 'o=s' =>
    \$outfile) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);



__END__

=head1

Normalizes matrix

=head1 USAGE

./code.pl [--i INPUT_FILE] [--p OPTION] [--o OUTPUT_FILE] [--helpp]

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

=item C<--helpp>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2011 May 15

=cut

