#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::Distributions;
use Time::SoFar qw(runtime);

sub splitgmt;

my $hubber = '/Genomics/ogtr04/arjunk/bin/Hubber';

my ($help, @idab, $infilt, $igmt, $igfilt, $igsfilt, $ipairfilt, $otab, $time, @p);
my $iming = 5; my $imaxg = 200;
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions('help' => \$help,
      'idab=s{,}' => \@idab,
       'infilt=s' => \$infilt,
         'igmt=s' => \$igmt,
       'igfilt=s' => \$igfilt,
        'iming=i' => \$iming,
        'imaxg=i' => \$imaxg,
      'igsfilt=s' => \$igsfilt,
    'ipairfilt=s' => \$ipairfilt,
         'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

my %sel_net = ();
if($infilt) {
    open DF, "$infilt" or die "Can't open $infilt!";
    while (<DF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $sel_net{$p[0]}++; }
    close DF; }

my %sel_pairs = ();
if($ipairfilt) {
    open PF, "$ipairfilt" or die "Can't open $ipairfilt!";
    while (<PF>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $p[2] = './contexts/'.$p[2];
        push(@{$sel_pairs{$p[1]}}, $p[2]); }
    close PF; }


open HH, ">$otab";
print HH "net\tid\tterm\t\tsize\tnum.in\tsum.in\tavg.in\t";
print HH "num.out\tsum.out\tavg.out\t";
print HH "num.inout\tsum.inout\tavg.inout\t";
print HH "segrn\tavg.segrn\n";


my %gs_pid;
if($igsfilt) {
    if($igfilt) {
        %gs_pid = splitgmt($igmt, $igsfilt, $igfilt); }
    else {
        %gs_pid = splitgmt($igmt, $igsfilt); } }
else {
    if($igfilt) {
        %gs_pid = splitgmt($igmt, $igfilt); }
    else {
        %gs_pid = splitgmt($igmt); } }
 

# my @gs_array = keys %$gs_desc;

my ($num_in, $num_inout, $sum_in, $sum_inout);
my ($gs, $avg_in, $sd_in, $avg_inout, $segrn, $avg_segrn);
my %gs_net_seg = (); my $count = 0;

print "\n";
foreach my $dab (@idab) {
    (my $net = $dab) =~ s/\.dab$//g; $net =~ s/\.qdab$//g; $net =~ s/^.*\///g;
    if($infilt) { unless(exists $sel_net{$net}) { next; } }
    print "$net ...\n";
    
    (my $ohub = $dab) =~ s/\.q*dab/\.hubber\.out/g; $ohub =~ s/^.*\///g;
    $count ++; $time = runtime(); print "$time:\t$count\t$ohub\n";

    if($ipairfilt) {
        my $ctxts = join ' ', @{$sel_pairs{$net}};
        print "\t$ctxts\n";
        `Hubber -e $igfilt -i $dab $ctxts > $ohub`; }
    else {
        `Hubber -e $igfilt -i $dab ./contexts/* > $ohub`; }

    open FH, "$ohub" or die "Can't open $ohub!";
    chomp(my @f=<FH>); close FH;

    shift(@f); shift(@f);
    foreach (@f) {
        @p = split '\t', $_;
        $gs = $gs_pid{$p[0]};

        $num_in = $p[7];
        $sd_in = $p[6];
        $avg_in = $p[5];
        $sum_in = ($avg_in*$num_in);

        $num_inout = ($p[4] - $num_in);
        $sum_inout = (($p[4]*$p[2]) - $sum_in);
        if($num_inout != 0) { $avg_inout = ($sum_inout/$num_inout); }
        else { $avg_inout = 'NA'; }

        $num_out = $num_inout - $num_in;
        $sd_out = $p[3];
        $sum_out = $sum_inout - $sum_in;
        if($num_out != 0) { $avg_out = ($sum_out/$num_out); }
        else { $avg_out = 'NA'; }

        if($sum_inout != 0) { $segrn = ($sum_in/$sum_inout); }
        else { $segrn = 'NA'; }

        if($avg_inout != 0) { $avg_segrn = ($avg_in/$avg_inout); }
        else { $avg_segrn = 'NA'; }

        print HH "$net\t$gs\t$p[0]\t$p[1]";
        print HH "\t$num_in\t", sprintf("%.3f\t%.3f", $sum_in, $avg_in);
        print HH "\t$num_out\t", sprintf("%.3f\t%.3f", $sum_out, $avg_out), "\t$num_inout";
        print HH sprintf("\t%.3f\t%.3f\t%.3f\t%.3f\n", $sum_inout, $avg_inout, $segrn, $avg_segrn);
    }
}
close HH;

$time = runtime(); print "\n$time: DONE\n\n";


# Subroutines
# ###########
sub splitgmt {
    my $gmt_file = shift;

    `mkdir -p ./contexts`;

    my %sel_gs = ();
    if($igsfilt) {
        my $gmt_filt = shift;
        open GS, "$gmt_filt" or die "Can't open $gmt_filt!";
        while (<GS>) {
            if($_ =~ /^#/) { next; }
            chomp($_); my @q = split '\t', $_;
            $sel_gs{$q[0]}++; }
        close GS; }

    my %sel_genes = ();
    if($igfilt) {
        my $gene_filt = shift;
        open GN, "$gene_filt" or die "Can't open $gene_filt!";
        while (<GN>) {
            if($_ =~ /^#/) { next; }
            chomp($_); my @q = split '\t', $_;
            $sel_genes{$q[0]}++; }
        close GN; }

    # my %gs_desc = my %gs_pid = ();
    my %gs_pid = ();
    
    open GMT, "$gmt_file" or die "Can't open $gmt_file!";
    while(<GMT>) {
        chomp($_); my @q = split '\t', $_;
        if($igsfilt) { unless(exists $sel_gs{$q[0]}) { next; } }

        my @genes = ();
        if($igfilt) {
            @genes = grep {exists $sel_genes{$_}} @q[2..$#q]; }
        else { @genes = @q[2..$#q]; }
        my $ngenes = scalar @genes;
        if(($ngenes < $iming) or ($ngenes > $imaxg)) { next; }

        (my $pid = $q[1]) =~ s/ \([0-9]*\)$//g;
        $pid = s/[ \/]/_/g; $pid =~ s/[;,'".]//g;
        $gs_pid{$pid} = $q[0]; shift(@q); shift(@q);

        open GS, ">./contexts/$pid";
        foreach my $g (@q) { print GS "$g\n"; }
        close GS; }
    close GMT;

    # return (\%gs_desc, \%gs_pid);
    return %gs_pid;
}


__END__

=head1

Get cohesiveness of all genesets in the network.

=head1 USAGE

run_hubber.pl [--idab NET_DAB] [--igmt GENESETS_GMT] [--infilt NET_FILTER]
[--gfilt GENE_FILTER] [--igsfilt GENESET_FILTER] [--ipairfilt NET-GENESET_FILTER] [--otab OUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes a network and a gmt file as inputs and calculates the
cohesiveness of genes in each geneset in the network. Cohesiveness is a measure
of how connected genes in a geneset are to each other relative to how they are
connected to all the genes in the network. Highly modular genesets will have
high cohesiveness.

=head1 ARGUMENTS

=over 12

=item C<--idab>

Network file in DAT/DAB format. Multiple DABs can be given using wildcard
matching.

=item C<--igmt>

Genesets file in GMT format.

=item C<--infilt>

(Optional) List of networks to consider for analysis.

=item C<--igfilt>

(Optional) List of genes to filter the genesets and network.

=item C<--igsfilt>

(Optional) List of genesets to consider for analysis.

=item C<--ipairfilt>

(Optional) List of network-geneset pairs to consider for analysis.

=item C<--otab>

Output table with each geneset per line.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 25

=cut

