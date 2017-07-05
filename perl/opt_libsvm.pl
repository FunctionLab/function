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
my $train = '/home/arjunk/software/svm/libsvm-3.12/svm-train';

# User Input
my(@issv, $help, $iker, $icost); my $icvk = 3; my $icvn = 5; my $imes = 'wauprc';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<1);
GetOptions( 'help' => \$help,
       'issv=s{,}' => \@issv,
          'iker=i' => \$iker,
          'icvk=i' => \$icvk,
          'icvn=i' => \$icvn,
          'imes=s' => \$imes,
         'icost=f' => \$icost) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

$iker = 2; $icvn = 1;

# For each dataset, for each paramter: n-times k-fold CV
my ($tmpf, @p, $score, $opt_c, $out_par, $time);
# my ($mean, $tempmean, $sqsum, $tempsqs, $ratio, $tempratio);
my (%nscore, $median, $mad, $tempratio, $ratio);

my @cost_array = ();
if($icost) { push(@cost_array, $icost); }
else {
    # for(my $log2c=-10; $log2c<=10; $log2c++) { # Could try [-10,10] or [-5,10]
    #     push(@cost_array, 2**$log2c); }
    for(my $log10c=-3; $log10c<=3; $log10c++) {
        push(@cost_array, 10**$log10c); }
}

my @gamma_array = ();
for(my $log10g=-5; $log10g<=1; $log10g++) {
    push(@gamma_array, 10**$log10g); }

my @cw_array = qw(0.01 0.05 0.1 0.25 0.5 0.75 1); # Class weights

my $m;
if($imes eq 'auc') { $m = 0; }
elsif($imes eq 'p10r') { $m = 1; }
elsif($imes eq 'p20r') { $m = 2; }
elsif($imes eq 'p50r') { $m = 3; }
elsif($imes eq 'bestf') { $m = 4; }
elsif($imes eq 'auprc') { $m = 5; }
elsif($imes eq 'wauprc') { $m = 6; }

# print "\nData\t#P\t#N\tC\tG\tW\tAUC\tP10R\tP20R\tP50R\tBest.F\tAUPRC\tWAUPRC\n";

foreach my $ssv (@issv) {
    $out_par = $ssv.'.par'; print "\n$ssv ...\n";
    open PAR, ">$out_par";
    print PAR "# $ssv\n#P\t#N\tC\tG\tW\tAUC\tP10R\tP20R\tP50R\tBest.F\tAUPRC\tWAUPRC\n";
    
    $opt_c = -1; $ratio = 0; $score = 0;
    if($icost) { $opt_c = $icost; }

    $tmpf = $out_par.'.tmp'; `rm -f $tmpf`;

    foreach my $c (@cost_array) {
        # print PAR sprintf("%.3g", $c);
        print "\tcost: ", sprintf("%.3g", $c);

        foreach my $g (@gamma_array) {
            # print PAR sprintf("\t%.3g", $g);
            print ", gamma: ", sprintf("%.3g", $g);

            foreach my $w (@cw_array) {
                # print PAR sprintf("\t%.3g", $w);
                print ", weight: ", sprintf("%.3g", $w), "\n\t\trep:";

                $median = 0; $mad = 0;

                %nscore = ();
                for(my $n=0; $n<$icvn; $n++) {
                    print " $n";
                    `$train -t $iker -c $c -g $g -v $icvk -e 0.001 -h 0 -w1 1 -w-1 $w $ssv | tail -9 > $tmpf`;

                    open TMP, "$tmpf"; chomp(my @tmp=<TMP>); close TMP;

                    @p = split ' ', $tmp[0]; print PAR "$p[2]";
                    @p = split ' ', $tmp[1]; print PAR "\t$p[2]\t$c\t$g\t$w";
                    foreach my $l (@tmp[2..8]) {
                        @p = split ' ', $l; print PAR "\t", sprintf("%.3f", $p[2]); }
                    print PAR "\n";

                    @p = split ' ', $tmp[$m];
                    $nscore{$n} = $p[2];

                    `rm -f $tmpf`;
                }

                $median = get_median(\%nscore);
                # print "\n\tscore: ", sprintf("%.2f", $median), "\n";
                # print PAR "\t", sprintf("%.2f", $median);
                if($icvn > 1) {
                    $mad = get_mad(\%nscore, $median);
                    $tempratio = ($median/$mad);
                    print PAR "\t", sprintf("%.2f", $mad);
                    print PAR "\t", sprintf("%.2f", $tempratio), "\n"; }

                # if($median > $score) {
                #     $score = $median; $opt_c = $c; }
            }
        }
    }

    close PAR;
    
    # ($tmpf = $ssv) =~ s/\.ssv/\.c$opt_c\.dck/g;
    # `$train -c $opt_c -v $icvk -e 0.001 $ssv > $tmpf`;
    print "\n\n"; }


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

=item C<--iker>

0 -- linear: u'*v

1 -- polynomial: (gamma*u'*v + coef0)^degree

2 -- radial basis function: exp(-gamma*|u-v|^2) (Default)

3 -- sigmoid: tanh(gamma*u'*v + coef0)

4 -- precomputed kernel (kernel values in training_set_file)

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

