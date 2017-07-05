#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);


my ($help, $imat, $imap, $ifiltc, $ifiltr, $itrimr, $ilogm, $omat);
my $imeasr = 'median';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions(     'help' => \$help,
              'imat=s' => \$imat,
              'imap=s' => \$imap,
              'itrimr' => \$itrimr,
            'imeasr=s' => \$imeasr,
               'ilogm' => \$ilogm,
            'ifiltr=s' => \$ifiltr,
            'ifiltc=s' => \$ifiltc,
              'omat=s' => \$omat        ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);


my %oid2nid = my %nid2oid = ();
open MAP, "$imap" or die "Can't open $imap!";
while (<MAP>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;
    if(($p[0] eq '') or ($p[0] =~ /\?/) or ($p[0] eq 'NA')) { next; }
    if(($p[1] eq '') or ($p[1] =~ /\?/) or ($p[1] eq 'NA')) { next; }

    my $oid = shift @p;
    if($itrimr) { $oid =~ s/\.[0-9][0-9]*$//g; }

    foreach my $nid (@p) {
        $oid2nid{$oid}{$nid}++;
        $nid2oid{$nid}{$oid}++; } }
close MAP;


my %filtc = ();
if($ifiltc) {
    open FC, "$ifiltc" or die "Can't open $ifiltc!";
    while (<FC>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $filtc{$p[0]}++; }
    close FC; }


my $tot_oid = my $num_oid_womap = 0;
my $num_oid_wmmap = my $num_oid_wumap = 0;
my $nline = 0; my %row_col_val = my @acol = ();

open OMAT, ">$omat";
open IMAT, "$imat" or die "Can't open $imat!";
while (<IMAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    $nline++;
    if($nline == 1) {
        shift @p; @acol = @p;
        #print OMAT join "\t", ('newid', @acol); print OMAT "\n";
        print OMAT "newid";
        foreach my $col (@acol) {
            if($ifiltc) {
                unless(exists $filtc{$col}) { next; }
                print OMAT "\t$col"; }
            else { print OMAT "\t$col"; } }
        print OMAT "\n";
        next; }

    my $oid = shift @p; $tot_oid++;
    if($itrimr) { $oid =~ s/\.[0-9][0-9]*$//g; }

    unless(exists $oid2nid{$oid}) {
        $num_oid_womap++; next; }
    if(scalar keys %{$oid2nid{$oid}} > 1) {
        $num_oid_wmmap++; next; }
    $num_oid_wumap++;

    for(my $j=0; $j<=$#p; $j++) {
        if($ifiltc) { unless(exists $filtc{$acol[$j]}) { next; } }
        unless(looks_like_number($p[$j])) { next; }
        foreach my $nid (keys %{$oid2nid{$oid}}) {
            push(@{$row_col_val{$nid}{$acol[$j]}}, $p[$j]); } } }
close IMAT;


my %filtr = ();
if($ifiltr) {
    open FR, "$ifiltr" or die "Can't open $ifiltr!";
    while (<FR>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $filtr{$p[0]}++; }
    close FR; }


my @anid = sort keys %row_col_val;
my $num_nid_wumap = my $num_nid_wmmap = 0;

ROW: foreach my $nid (@anid) {
    if($ifiltr) {
        unless(exists $filtr{$nid}) { next ROW; } }
    print OMAT "$nid";

    COL: foreach my $col (@acol) {
        if($ifiltc) {
            unless(exists $filtc{$col}) { next COL; } }

        if(exists $row_col_val{$nid}{$col}) {
            my $n = scalar @{$row_col_val{$nid}{$col}};
            my $avg_val;
            if($imeasr eq 'median') {
                $avg_val = median($row_col_val{$nid}{$col}); }
            elsif($imeasr eq 'mean') {
                $avg_val = mean($row_col_val{$nid}{$col}); }
            elsif($imeasr eq 'max') {
                $avg_val = max($row_col_val{$nid}{$col}); }
            elsif($imeasr eq 'sum') {
                $avg_val = sum($row_col_val{$nid}{$col}); }

            if($n == 1) { $num_nid_wumap++; }
            else { $num_nid_wmmap++; }

            #if($avg_val == 0) { print OMAT "\t0"; }
            #elsif($avg_val == 1) { print OMAT "\t1"; }
            #else {
            printf OMAT "\t%.5f", $avg_val; }
        else {
            print OMAT "\tNA"; } }
    print OMAT "\n"; }

print "\n> $imat\n\nNo. of cols: ", scalar @acol; print "\n\nTot. no. old ids: $tot_oid\n";
print "No. old ids w/ no mapping: $num_oid_womap\n";
print "No. old ids w/ many mappings: $num_oid_wmmap\n";
print "No. old ids w/ uniq mapping: $num_oid_wumap\n\n";

print "Tot. no. new ids: ", scalar @anid; print "\n";
print "No. new ids w/ many mappings: ", ($num_nid_wmmap/(scalar @acol)), "\n";
print "No. new ids w/ uniq mapping: ", ($num_nid_wumap/(scalar @acol)), "\n\n";

close OMAT;


# Calcualte mac
sub max {
    my $aref = shift;

    my @array = sort {$b <=> $a} grep {!/^NA$/} @$aref;
    return $array[0]; }


# Calculate mean
sub mean {
    my $aref = shift;
    
    if((scalar @$aref) == 1) {
        return ${$aref}[0]; }

    else {
        my $mean = my $n = 0;
        foreach (@$aref) {
            $n++;
            my $v = $_; if($ilogm) { $v = log($v + 0.0001); }
            $mean += ($v - $mean)/$n; }

        if($ilogm) { return (exp($mean) - 0.0001); }
        else { return $mean; } } }


# Calculates median
sub median {
    my $aref = shift;

    if((scalar @$aref) == 1) {
        return ${$aref}[0]; }

    else {
        my @array = sort {$a <=> $b} grep {!/^NA$/} @$aref;
        if(@array % 2) { return $array[@array/2]; }
        else { return (($array[(@array/2)-1] + $array[@array/2]) / 2); } } }


# Calculate sum
sub sum {
    my $aref = shift;

    my $sum = 0;
    foreach (@$aref) { $sum += $_; }

    return $sum; }


__END__

=head1

Map simple matrix from one row-identifier space to another.

=head1 USAGE

idmap_mat.pl [--imat INPUT_MATRIX] [--imap ID-MAP] [--ifiltr ROWS-TO-FILTER]
[--ifiltc COLUMNS-TO-FILTER] [--omat OUTPUT_MATRIX] [--help]

=head1 DESCRIPTION

This script will map a matrix from one row id-space to another. When there is
many-to-many mapping, all the old ids that map to multiple new ids will be
removed, and when multiple given ids map to a single new id, the rows of the old
ids will be averaged column-wise and assigned to the new id.

=head1 ARGUMENTS

=over 12

=item C<--imat>

Input matrix in a simple tab-separated format with header row containing column
identifiers, first columns containing row identifiers, and the rest of the file
filled with values.

=item C<--imap>

Mapping of row identifiers from one idspace to another. Tab-separated values of
<id1> <ide2>. Could contain many-to-many mapping, in one or many lines.

=item C<--ifiltr>

(Optional) List of row identifiers (in the new id space) that should be retained
in the output file.

=item C<--ifiltc>

(Optional) List of column identifiers that should be retained in the output file.

=item C<--itrimr>

(Optional) Trim row identifiers off numerical suffixes before mapping
(e.g. ENSG00234234.1 to ENSG00234234).

=item C<--ilogm>

(Optional) Calculate mean in log-space and unlog in the output.

=item C<--omat>

Output matrix in the new row-identifier space.

=item C<--help>

Prints this help message.

=back

=head1 PARAMETERS

=over 12

=item C<--p>

Parameter

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

