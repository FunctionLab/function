#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($help, $idat, $iself, $idir, $odat);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions( 'help' => \$help,
          'idat=s' => \$idat,
           'iself' => \$iself,
            'idir' => \$idir,
          'odat=s' => \$odat) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my (@p, $e);
my $ntote = 0; my $nselfe = 0;
my %uniq_genes = (); my %uniq_edges = ();

open FH, "$idat" or die "Can't open $idat!";
while (<FH>) {
    if($_ =~ /^#/) { next; }
    chomp($_); @p = split '\t', $_;

	$ntote++;
    
    # if($p[0] eq $p[1]) { unless ($iself) { next; } } # Removing self-edges

    if($idir) {
        $e = join '__', ($p[0], $p[1]); }
    else {
        $e = join '__', sort($p[0], $p[1]); }

    if(exists $uniq_edges{$e} and ($p[2] < $uniq_edges{$e})) { next; }
    $uniq_edges{$e} = $p[2]; # Recording unique edges
}
close FH;

open HH,">$odat";
my $nuniqe = 0;
foreach my $e (keys %uniq_edges) {
    @p = split '__', $e;

    if($p[0] eq $p[1]) {
        $nselfe++;
        unless($iself) { next; } }

	$nuniqe++;
	print HH "$p[0]\t$p[1]\t$uniq_edges{$e}\n";

    $uniq_genes{$p[0]}++; $uniq_genes{$p[1]}++;
}
close HH;

print "\nTot. no. edges: $ntote\nNo. self edges: $nselfe\nNo. uniq edges: $nuniqe";
print " (No. uniq genes: ", scalar keys %uniq_genes, ")\n\n";

__END__

=head1

Cleanup DAT by removing self-edges and retaining only unique-edges.

=head1 USAGE

cleanup_dat.pl [--idat INPUT_DAT] [--idir OPTION] [--odat OUTPUT_DAT] [--help]

=head1 DESCRIPTION

This script takes in a DAT and cleans it up by removing self-edges and retaining
only unique edges. If the --idir option is provided, then the network is treated
as directed and edge uniqueness is dertermined accordingly (A-B is different
from B-A).

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input network in DAT format.

=item C<--iself>

(Optional) If provided, self-edges are allowed. If not, they are removed (default).

=item C<--idir>

(Optional) If provided, the input network is treated as a directed graph.

=item C<--odat>

Output network in DAT format.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

