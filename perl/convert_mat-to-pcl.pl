#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);

my ($help, @imat, $imap, $opcl, @p, $time);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
       'imat=s{,}' => \@imat,
          'imap=s' => \$imap,
          'opcl=s' => \$opcl ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my %id_map = ();
if($imap) {
    open MAP, "$imap" or die "Can't open $imap!";
    while (<MAP>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        if(exists $id_map{$p[0]}) {
            $id_map{$p[0]} .= '|'.$p[1]; }
        else { $id_map{$p[0]} = $p[1]; }
    }
    close MAP;
}

if($#imat > 0) { undef $opcl; }

FILE: foreach my $amat (@imat) {
    if($amat =~ /\.pcl$/) {
        (my $omat = $amat) =~ s/\.pcl$/\.mat/g;
        open MAT, ">$omat";
        $time = runtime(); print "\n$time: $amat\t$omat";

        open PCL, "$amat" or die "Can't open $amat!";
        while (<PCL>) {
            if($_ =~ /^#/) { next; }
            if($_ =~ /^EWEIGHT/) { next; }

            chomp($_); @p = split '\t', $_;
            splice(@p, 1, 2);

            print MAT shift(@p);
            foreach my $val (@p) { print MAT "\t$val"; }
            print MAT "\n";
        }
        close PCL; close MAT;
        next FILE;
    }

    unless($imap) {
        %id_map = ();
        `rm -f tmp.rowids; cut -f1 $amat > tmp.rowids`;
        open MAP, "tmp.rowids"; chomp(my @ids = <MAP>); close MAP;
        shift(@ids); foreach (@ids) { $id_map{$_} = $_; }
    }

    open MAT, "$amat" or die "Can't open $amat!";
    chomp(my @mat = <MAT>); close MAT;

    unless($opcl) { ($opcl = $amat) =~ s/\.mat$/\.pcl/g; }
    open PCL, ">$opcl";
    $time = runtime(); print "\n$time: $amat\t$opcl";
    undef $opcl;

    my @header = split '\t', shift(@mat);
    print PCL "$header[0]\tName\tGWEIGHT";
    shift(@header); foreach (@header) { print PCL "\t$_"; }
    print PCL "\nEWEIGHT\t\t"; foreach (@header) { print PCL "\t1"; }
    print PCL "\n";

    my $newid;
    foreach (@mat) {
        @p = split '\t', $_;

        $newid = $p[0]; if(exists $id_map{$p[0]}) { $newid = $id_map{$p[0]}; }
        print PCL "$p[0]\t$newid\t1";

        shift(@p); # print "$newid\t", scalar @p, "\n";
        foreach (@p) {
            if($_ =~ /^$/) { $_ = 0; }
            print PCL "\t$_"; }
        print PCL "\n";
    }

    close PCL;
}

`rm -f tmp.rowids`;
$time = runtime(); print "\n$time: DONE\n\n";

__END__

=head1

Convert MAT to PCL (typically to use Sleipnir tools) or vice-versa.

=head1 USAGE

./convert_mat-to-pcl.pl [--imat INPUT_MAT/PCL] [--imap ID_MAPPING] [--opcl
OUTPUT_PCL] [--help]

=head1 DESCRIPTION

This script takes in a matrix in its simplest form - one header row with column
ids & the first column containing row ids - and converts it to PCL format. This
script can also be used to convert PCLs to MATs.

=head1 ARGUMENTS

=over 12

=item C<--imat>

Input matrix (simple format with one header row & row ids listed in column 1).
Multiple matrices can be provided using wildcard matching, if a mapping file is
also provided using the --imap option, then that mapping should apply to all the
matrices. Input can also be one or more PCL files, in which case, all other
options are ignored & the PCL are converted to MATs.

=item C<--imap>

(Optional) Gene/Row ID mapping file: <Given_ID> <New_ID>. Can have many-to-many
mapping.

=item C<--opcl>

Output matrix in PCL format.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2011 Oct 19

=cut

