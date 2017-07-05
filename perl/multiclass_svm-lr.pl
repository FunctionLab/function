#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);

sub get_cvcuts;
sub write_cvdat;
sub write_mc2bin;
sub train_test;
sub get_predscore;
sub print_pred;
sub dcheck;
sub sigmoid_train;
sub sigmoid_predict;

# Input
my(@idat, $imode, $ipar, $otab, $help); my ($icvn, $icvk) = (5, 3);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV<2);
GetOptions( 'help' => \$help,
       'idat=s{,}' => \@idat,
          'icvk=i' => \$icvk,
          'icvn=i' => \$icvn,
         'imode=s' => \$imode,
          'ipar=s' => \$ipar,
          'otab=s' => \$otab ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

# Output
open my $TAB, ">$otab";
print $TAB "#ID\tPAR\t#P\t#N\tAUC\tP10R\tP20R\tP50R\tAUPRC\n";

# Tools
my $subset = '/home/arjunk/software/svm/libsvm-3.11/tools/subset.py';
my $train_bin = '/home/arjunk/software/svm/liblinear-mod/train';
my $train_mul = '/home/arjunk/software/svm/liblinear-orig/train';
my $predict_bin = '/home/arjunk/software/svm/liblinear-mod/predict';
my $predict_mul = '/home/arjunk/software/svm/liblinear-orig/predict';

# Classification
my (@p, $nclass, %class_size, %dat_vec, $ngenes, %cv_idx, $y, $par);
my ($train_dat, $train_subdat, $test_dat, $test_subdat, $model, $A, $B);
my ($rclass, %orig_class, %tmp_score,
    %pred_score, $pred_out, $pred_labels, $pred_tab, $time);

my %dat_par = ();
if($ipar) {
    open PAR, "$ipar" or die "Can't open $ipar!";
    while (<PAR>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $dat_par{$p[0]} = $p[1];
    }
    close PAR;
}

# For each dataset ...
foreach my $ssv (@idat) {
    open DAT, "$ssv" or die "Can't open $ssv!";
        chomp(my @dat=<DAT>); close DAT;
    
    $ngenes = 0; %class_size = (); %dat_vec = (); # class -> cidx -> $vec
    if($ipar) { $par = $dat_par{$ssv}; } else { $par = 1; }
    
    # Hash: class -> cidx -> vec
    foreach my $vec (@dat) {
        $ngenes++;
        @p = split ' ', $vec;
        $class_size{$p[0]}++;
        $dat_vec{$p[0]}{$class_size{$p[0]}-1} = $vec;
    }

    $nclass = scalar keys %class_size;

    ($pred_labels = $ssv) =~ s/_[a-z]*\.[a-z]*$/\.pred.labels/g;
    `rm -f $pred_labels`;

    # For repeat run of CV ...
    for(my $n=0; $n<$icvn; $n++) {

        %cv_idx = (); # class -> fold -> @cidx
        foreach my $c (keys %class_size) {
            $cv_idx{$c} = get_cvcuts($dat_vec{$c}, $class_size{$c}, $icvk);
        }

        # For each CV fold ...
        for(my $k=0; $k<$icvk; $k++) {
            ($train_dat = $ssv) =~ s/_[a-z]*\.[a-z]*$/\.cv$k\.train/g;
            ($test_dat = $ssv) =~ s/_[a-z]*\.[a-z]*$/\.cv$k\.test/g;

            open my $TR, ">$train_dat";
            open my $TE, ">$test_dat";

            $y = 0; %orig_class = (); # prt_cidx -> $class
            foreach my $c (keys %class_size) {
                ($rclass, $y) =
                    write_cvdat($cv_idx{$c}, $dat_vec{$c}, $k, $TR, $TE, $y);
                @orig_class{keys %$rclass} = values %$rclass;
            }
            close $TR; close $TE;

            # MC - SVM / LR - One-vs-Rest
            if(($imode eq 'mc1') or ($imode eq 'mc2')) {
                ($pred_out = $test_dat) =~ s/$/\.pred/g;

                %pred_score = (); # prt_cidx -> class -> $score
                foreach my $c (keys %class_size) {
                    ($train_subdat = $train_dat) =~ s/\.train/\.sub-train\.c$c/g;
                    ($test_subdat = $test_dat) =~ s/\.test/\.sub-test\.c$c/g;

                    write_mc2bin($train_dat, $test_dat, $c);

                    if($imode eq 'mc1') {
                        ($A, $B) = train_test($train_subdat, $test_subdat, $c, $pred_out);
                        %tmp_score = get_predscore($pred_out, $A, $B); }

                    elsif($imode eq 'mc2') {
                        train_test($train_subdat, $test_subdat, $c, $pred_out);
                        %tmp_score = get_predscore($pred_out, $c); }

                    foreach my $idx (keys %tmp_score) {
                        $pred_score{$idx}{$c} = $tmp_score{$idx}; }
                }

                print_pred(\%pred_score, \%orig_class, $pred_labels);
            }

            # MC - SVM / LR - One-vs-One
            elsif(($imode eq 'mc3') or ($imode eq 'mc4')) {
                ($pred_out = $test_dat) =~ s/$/\.pred/g;

                %pred_score = (); # prt_cidx -> class -> $score

                my @class_array = sort keys %class_size; my ($ci, $cj);
                for(my $i=0; $i<$#class_array; $i++) {
                    $ci = $class_array[$i];

                    for(my $j=$i+1; $j<=$#class_array; $j++) {
                        $cj = $class_array[$j];

                        ($train_subdat = $train_dat)
                            =~ s/\.train/\.sub-train\.c$ci-c$cj/g;
                        ($test_subdat = $test_dat)
                            =~ s/\.test/\.sub-test\.c$ci-c$cj/g;

                        write_mc2bin($train_dat, $test_dat, $ci, $cj);

                        if($imode eq 'mc3') {
                            ($A, $B) = train_test($train_subdat, $test_subdat, $ci, $pred_out);
                            %tmp_score = get_predscore($pred_out, $A, $B); }
                        
                        elsif($imode eq 'mc4') {
                            train_test($train_subdat, $test_subdat, $ci, $pred_out);
                            %tmp_score = get_predscore($pred_out); }

                        foreach my $idx (sort {$a<=>$b} keys %tmp_score) {
                            $pred_score{$idx}{$ci} += $tmp_score{$idx};
                            $pred_score{$idx}{$cj} += (1-$tmp_score{$idx});
                        }
                    }
                }

                foreach my $idx (keys %pred_score) {
                    foreach my $c (@class_array) {
                        $pred_score{$idx}{$c} /= scalar @class_array;
                    }
                }

                print_pred(\%pred_score, \%orig_class, $pred_labels);
            }

            # MC - Crammer & Singer
            elsif($imode eq 'mc5') {
                ($model = $train_dat) =~ s/$/\.model/g;
                `$train_mul -c $par -e 0.001 -s 4 $train_dat $model`;
                `$predict_mul $test_dat $model tmp.pred`;
                `cat tmp.pred >> $pred_labels; rm -f tmp.pred $model`;
            }

            # BI - SVM / LR - One-vs-Rest
            elsif(($imode eq 'bi1') or ($imode eq 'bi2')) {
                my $c = 1; %pred_score = ();
                ($train_subdat = $train_dat) =~ s/\.train/\.sub-train\.c$c/g;
                ($pred_out = $test_dat) =~ s/$/\.pred/g;

                write_mc2bin($train_dat, $test_dat, $c);

                if($imode eq 'bi1') {
                    ($A, $B) = train_test($train_subdat, $test_dat, $c, $pred_out);
                    %tmp_score = get_predscore($pred_out, $A, $B); }

                elsif($imode eq 'bi2') {
                    train_test($train_subdat, $test_dat, $c, $pred_out);
                    %tmp_score = get_predscore($pred_out); }

                foreach my $idx (keys %tmp_score) {
                    $pred_score{$idx}{$c} = $tmp_score{$idx}; }

                print_pred(\%pred_score, \%orig_class, $pred_labels);
            }

            # BI - SVM / LR - One-vs-One
            elsif(($imode eq 'bi3') or ($imode eq 'bi4')) {
                ($pred_out = $test_dat) =~ s/$/\.pred/g;

                %pred_score = (); # prt_cidx -> class -> $score

                my @class_array = sort keys %class_size; my ($ci, $cj);
                $ci = 1; my $i = 0;

                for(my $j=$i+1; $j<=$#class_array; $j++) {
                    $cj = $class_array[$j];

                    ($train_subdat = $train_dat)
                        =~ s/\.train/\.sub-train\.c$ci-c$cj/g;
                    ($test_subdat = $test_dat)
                        =~ s/\.test/\.sub-test\.c$ci-c$cj/g;

                    write_mc2bin($train_dat, $test_dat, $ci, $cj);

                    if($imode eq 'bi3') {
                        ($A, $B) = train_test($train_subdat, $test_subdat, $ci, $pred_out);
                        %tmp_score = get_predscore($pred_out, $A, $B); }

                    elsif($imode eq 'bi4') {
                        train_test($train_subdat, $test_subdat, $ci, $pred_out);
                        %tmp_score = get_predscore($pred_out); }

                    foreach my $idx (sort {$a<=>$b} keys %tmp_score) {
                        $pred_score{$idx}{$ci} += $tmp_score{$idx};
                    }
                }

                foreach my $idx (keys %pred_score) {
                    $pred_score{$idx}{$ci} /= scalar @class_array;
                }

                print_pred(\%pred_score, \%orig_class, $pred_labels);
            }

            `rm -f $train_dat $test_dat`;
        }
    }

    ($pred_tab = $ssv) =~ s/_[a-z]*\.[a-z]*$/\.$imode\.pred.txt/g;
    dcheck($pred_labels, $pred_tab, $par, $TAB);
    
    `rm -f *.sub-train.* *.sub-test.*`;

    $ssv =~ s/_[a-z]*\.[a-z]*$//g;
    $time = runtime(); print "\n$time: $ssv\t$ngenes";
}

print "\nDONE\n\n";

# Given indices and size for a class, assigns indices to cv-fold
sub get_cvcuts {
    my $dat = shift;
    my $size = shift;
    my $cvk = shift;

    my $cvcut = int($size/$cvk + 0.5);

    # Randomize the indices
    my @rnd_idx = shuffle(keys %$dat);
    my %cvi = ();
    my $i = 0; my $s = 0;

    # For each fold ...
    for(my $k=0; $k<$cvk; $k++) {
        $s += $cvcut;
        if($k == ($cvk-1)) { $s = $size; }

        # Record the random index for that fold
        do{ push(@{$cvi{$k}}, $rnd_idx[$i]);
            $i++;
        } while ($i < $s);
    }

    return \%cvi;
}

# Given indices and data for a class, prints training & testing data
sub write_cvdat {
    my $idx = shift;
    my $dat = shift;
    my $k = shift;
    my $tr = shift;
    my $te = shift;
    my $x = shift;

    my @q; my %oc = ();

    # Print examples split into train and test
    foreach my $f (keys %$idx) {
        if($f == $k) {
            foreach my $i (@{$idx->{$f}}) {
                print $te $dat->{$i}, "\n";
                @q = split ' ', $dat->{$i};
                $oc{$x} = $q[0]; $x++;
            }
        }
        else {
            foreach my $i (@{$idx->{$f}}) {
                print $tr $dat->{$i}, "\n";
            }
        }
    }

    return (\%oc, $x);
}

# Given training data & a class, converts labels of examples with that class to
# 1 & the others to -1
sub write_mc2bin {
    my $trdat = shift;
    my $tedat = shift;

    my $subtr = $trdat;
    my $subte = $tedat;
    my $vec;

    if(($imode eq 'mc1') or ($imode eq 'mc2')
        or ($imode eq 'bi1') or ($imode eq 'bi2')) {
        my $class = shift;
        $subtr =~ s/\.train/\.sub-train\.c$class/g;
        $subte =~ s/\.test/\.sub-test\.c$class/g;

        open TR, "$trdat"; open STR, ">$subtr";
        while(<TR>) {
            chomp($vec = $_);
            if($vec =~ /^$class /) { $vec =~ s/^$class /+1 /g; }
            else { $vec =~ s/^[0-9] /-1 /g; }
            print STR "$vec\n";
        }
        close TR; close STR;

        open TE, "$tedat"; open STE, ">$subte";
        while(<TE>) {
            chomp($vec = $_);
            if($vec =~ /^$class /) { $vec =~ s/^$class /+1 /g; }
            else { $vec =~ s/^[0-9] /-1 /g; }
            print STE "$vec\n";
        }
        close TE; close STE;
    }

    elsif(($imode eq 'mc3') or ($imode eq 'mc4')
        or ($imode eq 'bi3') or ($imode eq 'bi4')) {
        my $ci = shift;
        my $cj = shift;
        $subtr =~ s/\.train/\.sub-train\.c$ci-c$cj/g;
        $subte =~ s/\.test/\.sub-test\.c$ci-c$cj/g;

        open TR, "$trdat"; open STR, ">$subtr";
        while(<TR>) {
            chomp($vec = $_);
            if($vec =~ /^$ci /) {
                $vec =~ s/^$ci /+1 /g; print STR "$vec\n"; }
            elsif($vec =~ /^$cj /) {
                $vec =~ s/^$cj /-1 /g; print STR "$vec\n"; } 
        }
        close TR; close STR;

        open TE, "$tedat"; open STE, ">$subte";
        while(<TE>) {
            chomp($vec = $_);
            if($vec =~ /^$ci /) {
                $vec =~ s/^$ci /+1 /g; print STE "$vec\n"; }
            elsif($vec =~ /^$cj /) {
                $vec =~ s/^$cj /-1 /g; print STE "$vec\n"; }
            else {
                $vec =~ s/^[0-9]* /0 /g; print STE "$vec\n"; }
        }
        close TE; close STE;
    }
}

# Given training, testing data, runs training, testing & writes out the
# predictions
sub train_test {
    my $trdat = shift;
    my $tedat = shift;
    my $class = shift;
    my $prdat = shift;

    my $model = $trdat.'.model';
    my @q; my %oclass = (); my %pscore = ();

    if(($imode eq 'mc1') or ($imode eq 'bi1')
        or ($imode eq 'mc3') or ($imode eq 'bi3')) {
        `$train_bin -c $par -e 0.001 $trdat`;
        `$predict_bin $trdat $model tmp.pred > $prdat`;

        open PR, "$prdat";
        while(<PR>) {
            if($_ =~ /^#/) { next; }
            if($_ =~ /^labels/) { next; }
            if($_ =~ /^Accuracy/) { next; }
            chomp($_); @q = split '\t', $_;

            if($q[1] eq $class) { $oclass{$q[0]} = 1; }
            else { $oclass{$q[0]} = -1; }

            $pscore{$q[0]} = $q[2];
        }
        my ($A, $B) = sigmoid_train(\%pscore, \%oclass);

        close PR; `rm -f $prdat`;

        `$predict_bin $tedat $model tmp.pred > $prdat`;
        `rm -f $trdat $model tmp.pred`;

        return ($A, $B);
    }

    elsif(($imode eq 'mc2') or ($imode eq 'bi2')
        or ($imode eq 'mc4') or ($imode eq 'bi4')) {
        `$train_mul -c $par -e 0.001 -s 0 $trdat`;
        `$predict_mul -b 1 $tedat $model $prdat`;
        `rm -f $trdat $model`;
    }
}

# Given prediction data and a class, estimates prediction score/probability for
# all the predicted examples for that class
sub get_predscore {
    my $prdat = shift;
    my ($A, $B) = @_;

    my %pscore = ();
    my (@q, $more, $less);
    
    open PR, "$prdat";
    while(<PR>) {
        if($_ =~ /^#/) { next; }
        if($_ =~ /^labels/) { next; }
        chomp($_); @q = split '\t', $_;

        if(($imode eq 'mc1') or ($imode eq 'bi1')
            or ($imode eq 'mc3') or ($imode eq 'bi3')) {
            $pscore{$q[0]} = sigmoid_predict($q[2], $A, $B); }

        elsif(($imode eq 'mc2') or ($imode eq 'bi2')
            or ($imode eq 'mc4') or ($imode eq 'bi4')) {
            ($less, $more) = sort {$a<=>$b} ($q[3], $q[4]);
            if($q[2] == 1) { $pscore{$q[0]} = $more; }
            elsif($q[2] == -1) { $pscore{$q[0]} = $less; }
        }
    }
    close PR;

    `rm -f $prdat`;
    return %pscore;
}

# Given prediction scores and original classes, writes out a table of original
# class, predicted class & the associated score
sub print_pred {
    my $pscore = shift;
    my $oc = shift;
    my $outlab = shift;

    my ($c, $sum, $max234, $ps);
    open PL, ">>$outlab";
    foreach my $idx (keys %$pscore) {
        if(($imode eq 'mc1') or ($imode eq 'mc2')
            or ($imode eq 'mc3') or ($imode eq 'mc4')) {
            $ps = -1000; $sum = 0.001; $max234 = 0;
            foreach $c (keys %{$pscore->{$idx}}) {
                $sum += $pscore->{$idx}{$c};
                if($c != 1) {
                    if($pscore->{$idx}{$c} > $max234) {
                        $max234 = $pscore->{$idx}{$c}; } }
            }
            $c = 1; $ps = $pscore->{$idx}{$c};
            # $ps = (-1)*$ps*log(1-($ps/$sum))/log(2);
            $ps = $ps*log($ps/$max234)/log(2);
        }
        elsif(($imode eq 'bi1') or ($imode eq 'bi2')
            or ($imode eq 'bi3') or ($imode eq 'bi4')) {
            $c = 1; $ps = $pscore->{$idx}{$c};
        }
        # print PL $oc->{$idx}, "\t$c\t", sprintf("%.3f", $ps), "\n";
        print PL $oc->{$idx}, "\t", sprintf("%.6f", $ps), "\n";
    }
    close PL;
}

sub dcheck {
    my $predl = shift;
    my $predt = shift;
    my $ppar = shift;
    my $FH = shift;

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

# Learns the sigmoid model from SVM decision values
# input: decision_values, real_labels{1,-1}, #pos_instances, #neg_instances
# output: [A, B] that minimize sigmoid likelihood
sub sigmoid_train {
    my $deci = shift;
    my $label = shift;
    my ($prior1, $prior0) = @_;

	#Count prior0 and prior1 if needed
    my $num_eg = scalar keys %$label;
    if((not defined $prior1) or (not defined $prior0)) {
        ($prior1, $prior0) = (0, 0);
        foreach my $i (keys %$label) {
            if($label->{$i} > 0) { $prior1++; }
            else { $prior0++; }
        }
    }

    #Parameter Setting
    my $max_iter = 100; # Max. no. iterations
    my $min_step = 1e-10; # Min. step taken in line search
    my $sigma = 1e-12; # For numerically strict PD of Hessian
    my $eps = 1e-5;

    #Construct Target Support
    my $hiTarget = ($prior1+1.0)/($prior0+2.0);
    my $loTarget = 1/($prior0+2.0);
    my %tar = ();

    foreach my $i (keys %$label) {
        if($label->{$i} > 0) { $tar{$i} = $hiTarget; }
        else { $tar{$i} = $loTarget; }
    }

    #Initial Point & Initial Fun Value
    my ($A, $B) = (0.0, log(($prior0+1.0)/($prior1+1.0)));
    my $fval = 0.0; my $fApB;

    foreach my $i (keys %$label) {
        $fApB = ($deci->{$i})*$A + $B;
        if($fApB >= 0) {
            $fval += $tar{$i}*$fApB + log(1+exp(-$fApB)); }
        else { $fval += ($tar{$i}-1)*$fApB + log(1+exp($fApB)); }
    }

    my ($itr, $h11, $h22, $h21, $g1, $g2, $p, $q, $d1, $d2);
    my ($det, $dA, $dB, $gd, $step_size, $newA, $newB, $newf);
    ITER: for($itr=0; $itr<$max_iter; $itr++) {
        #Update Gradient & Hessian (use H' = H + sigma I)
        $h11 = $h22 = $sigma; #Numerically ensures stricy PD
        $h21 = $g1 = $g2 = 0.0;
        foreach my $i (keys %$label) {
            $fApB = $deci->{$i}*$A + $B;
            if($fApB >= 0) {
                $p = exp(-$fApB)/(1.0+exp(-$fApB));
                $q = 1.0/(1.0+exp(-$fApB));
            }
            else {
                $p = 1.0/(1.0+exp($fApB));
                $q = exp($fApB)/(1.0+exp($fApB));
            }
            $d2 = $p*$q;
            $h11 += $deci->{$i}*$deci->{$i}*$d2;
            $h22 += $d2;
            $h21 += $deci->{$i}*$d2;
            $d1 = $tar{$i}-$p;
            $g1 += $deci->{$i}*$d1;
            $g2 += $d1;
        }

        #Stopping Criteria
        if((abs($g1) < $eps) and (abs($g2) < $eps)) {
            last ITER; }

        #Finding Newton direction: -inv(H') * g
        $det = $h11*$h22 - $h21*$h21;
        $dA = (-1)*($h22*$g1 - $h21*$g2) / $det;
        $dB = (-1)*(-$h21*$g1 + $h11*$g2) / $det;
        $gd = $g1*$dA + $g2*$dB;

        #Line Search
        $step_size = 1;
        STEP: while ($step_size >= $min_step) {
            $newA = $A + $step_size*$dA;
            $newB = $B + $step_size*$dB;

            #New function value
            $newf = 0.0;
            foreach my $i (keys %$label) {
                $fApB = $deci->{$i}*$newA + $newB;
                if($fApB >= 0) {
                    $newf += $tar{$i}*$fApB + log(1+exp(-$fApB)); }
                else { $newf += ($tar{$i}-1)*$fApB + log(1+exp($fApB)); }
            }

            #Check sufficient decrease
            if($newf < $fval + 0.0001*$step_size*$gd) {
                ($A, $B, $fval) = ($newA, $newB, $newf);
                last STEP;
            }
            else { $step_size = $step_size/2.0; }
        }

        if($step_size < $min_step) {
            print "line search fails: $A, $B, $g1, $g2, $dA, $dB, $gd\n";
            return ($A, $B);
        }
    }

    if($itr >= $max_iter) {
        print "reaching max. iterations: $g1, $g2"; }

    return ($A, $B);
}

# input: decision_values & Platt parameters [A, B]
# output: predicted probability
sub sigmoid_predict {
    my $deci = shift;
    my ($A, $B) = @_;

    my $fApB = $deci*$A + $B;
    if($fApB >= 0) {
        return exp(-$fApB)/(1.0+exp(-$fApB));  }
    else { return 1.0/(1+exp($fApB)); }
}

__END__

=head1

Cross-validation using LIBLINEAR for multi-class classification.

=head1 USAGE

cv_liblinear.pl [--idat INPUT_DATA_SSV] [--icvk CV_FOLD] [--icvn CV_NUM]
    [--imode LEARNING_MODE] [--ipar PARAMTER] [--otab OUTPUT_TABLE] [--help]

=head1 DESCRIPTION

This script takes ssv-format data with labels and features, and runs n-time
k-fold CV using LIBLINEAR to optimize the cost parameter and report prediction
performance (say, F1score). For binary classification, just use the opt-c code.
This one is designed to work for multi-class classification.

The most important purpose of this script is to run multi-class SVM, keeping in
mind that Class1 is most important. So, all optimizations & evaluations are done
based on Class1, i.e., [TP,TN,FP,FN] counts are calculated based on Class 1 vs.
Rest.

=head1 ARGUMENTS

=over 12

=item C<--idat>

Input datafiles in SSV format. Can give multiple files with wildcard matching.

=item C<--icvk>

No. of folds of CV

=item C<--icvn>

No. of times to run k-fold CV

=item C<--imode>

# Multi-class

[mc1] Multi-class SVM - Linear Kernel - One-vs-Rest (using SVM decision scores) (Winner-Takes-All)

[mc2] Multi-class Logistic Regression - One-vs-Rest

[mc3] Multi-class SVM - One-vs-One (Voting/Max-Wins)

[mc4] Multi-class Logistic Regression - One-vs-One

[mc5] Multi-class SVM - Direct method (Crammer & Singer)

[mc6] Multi-class SVM - DAGSVM

[mc7] Multi-class SVM - One-vs-One (Pairwise Coulping)

[mc8] Multi-class SVM - RBF Kernel - One-vs-Rest (using Platt's post. prob.) (Winner-Takes-All)

# Binary - Class1 vs. Rest - w/ stratified splitting of data (maintaining class distributions) for CV.

[bi1] Binary SVM - Linear Kernel - One-vs-Rest

[bi2] Binary Logistic Regression - One-vs-Rest

[bi3] Binary SVM - Linear Kernel - One-vs-One

[bi4] Binary Logistic Regression - One-vs-One

[bi5] Binary SVM - RBF Kernel - One-vs-Rest

=item C<--ipar>

Parameter for the chosen method. If the method SVM, then the parameter is the
cost. If this option is not provided, an optimal parameter is chosen using
cross-validation.

=item C<--otab>

Output table reporting the cost parameter and the best performance score.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 April 5

=cut

