#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;
use Time::SoFar qw(runtime);

sub get_median;
sub get_mad;

# Tools
my $train = '/home/arjunk/software/svm/liblinear-mod/train';

# User Input
my(@issv, $help, $icost); my $icvk = 3; my $icvn = 10; my $imes = 'wauprc';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<1);
GetOptions( 'help' => \$help,
       'issv=s{,}' => \@issv,
          'icvk=i' => \$icvk,
          'icvn=i' => \$icvn,
          'imes=s' => \$imes,
         'icost=f' => \$icost) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

# For each dataset, for each paramter: n-times k-fold CV
my ($tempfile, @p, $score, $opt_c, $out_par, $time);
# my ($mean, $tempmean, $sqsum, $tempsqs, $ratio, $tempratio);
my (%nscore, $median, $mad, $tempratio, $ratio);

my @cost_array = ();
if($icost) { push(@cost_array, $icost); }
else {
    # for(my $log2c=-10; $log2c<=10; $log2c++) { # Could try [-10,10] or [-5,10]
    #     push(@cost_array, 2**$log2c); }
    # for(my $log10c=-3; $log10c<=2; $log10c++) {
    #     push(@cost_array, 10**$log10c); }
    @cost_array = qw(0.001 0.01 0.1 1 10 100 1000);
}

my $m;
if($imes eq 'auc') { $m = 0; }
elsif($imes eq 'p10r') { $m = 1; }
elsif($imes eq 'p20r') { $m = 2; }
elsif($imes eq 'p50r') { $m = 3; }
elsif($imes eq 'bestf') { $m = 4; }
elsif($imes eq 'auprc') { $m = 5; }
elsif($imes eq 'wauprc') { $m = 6; }

# print "\nData\tOpt.C\t#P\t#N\tAUC\tP10R\tP20R\tP50R\tBest.F\tAUPRC\tWAUPRC\n";

foreach my $ssv (@issv) {
    ($out_par = $ssv) =~ s/ssv$/par/g; print "\n$ssv ...\n";
    open PAR, ">$out_par";
    
    $opt_c = -1; $ratio = 0; $score = 0;
    if($icost) { $opt_c = $icost; }

    foreach my $c (@cost_array) {
        # $mean = 0; $tempmean = 0;
        # $sqsum = 0; $tempsqs = 0;
        $median = 0; $mad = 0;
        $tempfile = $out_par.'.c';
        $time = runtime(); print "\t$time: cost: ", sprintf("%.3g", $c), "\trep:";
        print PAR sprintf("%.3g", $c);

        %nscore = ();
        for(my $n=0; $n<$icvn; $n++) {
            print " $n";
            `rm -f $tempfile`;
            # `$train -c $c -v $icvk -e 0.001 $ssv | tail -7 > $tempfile`;
            `$train -s 2 -c $c -v $icvk $ssv | tail -7 > $tempfile`;

            open TMP, "$tempfile"; chomp(my @tmp=<TMP>); close TMP;

            @p = split ' ', $tmp[$m];
            $nscore{$n} = $p[2];

            # $tempmean = $mean;
            # $mean = $tempmean + ($p[2]-$tempmean)/($n+1);
            # $tempsqs = $sqsum;
            # $sqsum = $tempsqs + ($p[2]-$tempmean)*($p[2]-$mean);
        }

        # $tempratio = ($mean/sqrt($sqsum/($icvn-1)));
        # print PAR "\tM: ", sprintf("%.2f", $mean);
        # print PAR "\tS: ", sprintf("%.2f", sqrt($sqsum/($icvn-1)));
        # print PAR "\tR: ", sprintf("%.2f", $tempratio), "\n";
        $median = get_median(\%nscore);
        $time = runtime(); print "\t$time: score: ", sprintf("%.2f", $median), "\n";
        print PAR "\t", sprintf("%.2f", $median);
        if($icvn > 1) {
            $mad = get_mad(\%nscore, $median);
            $tempratio = ($median/$mad);
            print PAR "\t", sprintf("%.2f", $mad);
            print PAR "\t", sprintf("%.2f", $tempratio); }
        print PAR "\n";

        # if($tempratio > $ratio) {
        #    $ratio = $tempratio; $score = $median; $opt_c = $c; }
        if($median > $score) {
            $score = $median; $opt_c = $c; }
    }

    `rm -f $tempfile`;
    # $time = runtime(); print "$time: $ssv: opt_c: $opt_c $p[0]: $score\n";
    # print "$ssv: opt_c: $opt_c $p[0]: $score\n";
    close PAR;
    
    ($tempfile = $ssv) =~ s/\.ssv/\.c$opt_c\.dck/g;
    # `$train -c $opt_c -v $icvk -e 0.001 $ssv | tail -9 > $tempfile`;
    # `$train -c $opt_c -v $icvk -e 0.001 $ssv > $tempfile`;
    `$train -s 2 -c $opt_c -v $icvk $ssv > $tempfile`;
    open TMP, "$tempfile"; chomp(my @tmp=<TMP>); close TMP;
    print "\nData\tOpt.C\t#P\t#N\tAUC\tP10R\tP20R\tP50R\tBest.F\tAUPRC\tWAUPRC\n";
    print "$ssv\t$opt_c";
    foreach my $l (@tmp[($#tmp-8)..($#tmp-7)]) {
        @p = split ' ', $l; print "\t$p[2]"; }
    foreach my $l (@tmp[($#tmp-6)..$#tmp]) {
        @p = split ' ', $l; print "\t", sprintf("%.3f", $p[2]); }
    print "\n\n"; }

$time = runtime(); print "$time: DONE\n\n";


sub get_median {
    my $href = shift;
    my @ary = sort {$a <=> $b} grep {!/^NA$/} (values %{$href});
    if(@ary % 2) { return $ary[@ary/2]; }
    else { return (($ary[(@ary/2)-1] + $ary[@ary/2]) / 2); }
}

sub get_mad {
    my $href = shift;
    my $med = shift;
    foreach (keys %{$href}) { $href->{$_} = abs($href->{$_} - $med); }
    my $mad = get_median($href);
    return $mad;
}

__END__

=head1 NAME

opt-c_liblinear.pl - Runs n-times k-fold CV and finds optimal cost parameter.

=head1 USAGE

opt-c_liblinear.pl [--issv INPUT_DATA_SSV] [--icvk CV_FOLD] [--icvn CV_NUM]
    [--help]

=head1 DESCRIPTION

This script takes ssv-format data with labels and features, and runs n-time
k-fold CV using LIBLINEAR to optimize the cost parameter and report prediction
performance (say, F1score).

=head1 ARGUMENTS

=over 12

=item C<--issv>

Input datafiles in SSV format. Can give multiple files with wildcard matching.

=item C<--icost>

Cost parameter. (Optional) If provided, only this value will be used. If not
provided, a range of values will be tried.

=item C<--icvk>

(Optional) No. of folds of CV. (Default = 3)

=item C<--icvn>

(Optional) No. of times to run k-fold CV. (Default = 10)

=item C<--imes>

(Optional) Performance measure to use to choose optimal cost: 'auc', 'p10r',
'p20r', 'p50r', 'bestf', 'auprc', or 'wauprc.'

=item C<--help>

Prints this documentation

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 March 5

=cut

