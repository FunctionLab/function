#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Pod::Usage;
#use Data::Dumper;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);

sub print_pos {
    my $genes_hashref = shift;
    my $fh1 = shift;
    # my $fh2 = shift;
    my %pos_genes = ();

    foreach my $g (keys %{$genes_hashref}) {
        print $fh1 "$g\n"; # print $fh2 "$g\n";
        $pos_genes{$g}++;
    }

    return \%pos_genes;
}

sub print_neg {
    my $npos = shift;
    my $rand_arrayref = shift;
    my $pos_hashref = shift;
    my $class_hashref = shift;
    my $fh = shift;
    my $nneg = 0;

    while($nneg < $npos) {
        my $g = shift(@{$rand_arrayref});
        if(exists $pos_hashref->{$g}) { next; }
        unless(exists $class_hashref->{$g}) {
            print $fh "$g\n";
            $nneg++;
        }
    }
}

sub write_svmdat {
    my $dat2dab = '/home/arjunk/src/c++/sleipnir/tools/Dat2Dab/Dat2Dab';

    my $net = shift;
    my $list1 = shift; my $list2 = shift;
    my $tag = shift;
    my $pos_hashref = shift;
    
    my $mat = $list1; $mat =~ s/_genes.txt/_$tag.mat/g;
    my $dat = $mat; $dat =~ s/\.mat/\.ssv/g;
    my (@q, $tempfile);

    `$dat2dab -i $net -t $list1 -T $list2 > $mat`;

    open MAT, "$mat"; chomp(my @dmat = <MAT>); close MAT;
    open SSV, ">$dat";

    shift(@dmat);
    foreach my $line (@dmat) {
        @q = split '\t', $line;
        
        if(exists $pos_hashref->{$q[0]}) { print SSV "+1"; }
        else { print SSV "-1"; }

        for(my $j=1; $j<=$#q; $j++) {
            print SSV " $j:$q[$j]";
        }
        print SSV "\n";
    }

    close SSV;
    `rm -f $mat`;

    `impute_missing_values.pl $dat`;
    $tempfile = $dat.'.imputed'; `mv $tempfile $dat`;

    return $dat;
}

sub svm_gridsearch {
    my $svm_train = '/home/arjunk/software/svm/libsvm-3.11/svm-train';
    my $grid_search = '/home/arjunk/software/svm/libsvm-3.11/tools/grid.py';
    my $gnuplot = '/home/arjunk/src/misc/gnuplot-4.4.4/src/gnuplot';
    my $grid_par = '-log2c -5,5,1 -log2g 1,1,1 -v 3 -svmtrain '.$svm_train.' -gnuplot '.$gnuplot.' -t 0';

    my $dat = shift;
    my $tempfile = $dat.'.grid';

    `python $grid_search $grid_par $dat > $tempfile`;

    open PH, "$tempfile"; chomp(my @svmpar=<PH>); close PH;
    my @q = split ' ', $svmpar[$#svmpar];

    return ($q[0], $q[2]);
}

sub svm_traintest {
    my $svm_train = '/home/arjunk/software/svm/libsvm-3.11/svm-train';
    my $svm_test = '/home/arjunk/software/svm/libsvm-3.11/svm-predict';
    
    my $train_dat = shift;
    my $cost = shift;
    my $test_dat = shift;

    `$svm_train -t 0 -c $cost $train_dat`;

    my $model = $train_dat.'.model';
    my $prediction = $test_dat.'.pred';
    my $tempfile = $test_dat.'.acc';

    `$svm_test $test_dat $model $prediction > $tempfile`;

    open AH, "$tempfile"; chomp(my $acc=<AH>); close AH;
    my @q = split ' ', $acc; $q[2] =~ s/%//g;;

    return $q[2];
}

# my $tissuenet_dir = '/Genomics/ogtr03/cgreene/tissue-networks/prop/pruned/standard_integration/predictions/';
my $tissuedat_dir = '/home/arjunk/data/gene-expression/cell-tissue-specificity/';
my $tissuenet_dir = '/bb/vol0/cgreene/tissue-networks/prop/pruned/average_proportions/predictions/';

# Tissue-specific genes *.gmt
# List of ubiquitous genes
# Disease gene association *.gmt
# Disease-tissue associations
# List of background genes in all tissue-networks
# Output: Table of selected diseases and their prediction accuracies

my $in_tspfc = $tissuedat_dir.'hprd/hprd-brenda/human_tissue-specific_gene-expression_pruned.gmt';
my $in_ubiq = $tissuedat_dir.'ubiquitous_genes/alldb_ubiquitous-genes_human.txt';
my $in_bglist = '/home/arjunk/data/networks/human/tissue_networks/genes_in_tissue-networks.txt';
my($in_dis, $in_dis_tis, $out_table) = @ARGV;

open TH, "$in_tspfc" or die "Can't open $in_tspfc!";
    chomp(my @tspfc=<TH>); close TH;
open UH, "$in_ubiq" or die "Can't open $in_ubiq!";
    chomp(my @ubiq=<UH>); close UH;
open DH, "$in_dis" or die "Can't open $in_dis!";
    chomp(my @dis=<DH>); close DH;
open DTH, "$in_dis_tis" or die "Can't open $in_dis_tis!";
    chomp(my @dt=<DTH>); close DTH;
open BG, "$in_bglist" or die "Can't open $in_bglist!";
    chomp(my @bg=<BG>); close BG;
open TB, ">$out_table";

my $time;
# Indexing background list of genes
my %bg_genes = ();
foreach (@bg) { $bg_genes{$_}++; }

# Indexing ubiquitous and tissue-specific genes
my %ubiq_gene = ();
foreach (@ubiq) { unless(exists $bg_genes{$_}) { next; } $ubiq_gene{$_}++; }

my %tissue_genes = (); my %tspfc_genes = (); my ($tissue, @p);
foreach (@tspfc) {
    @p = split '\t', $_;
    $tissue = shift(@p); shift(@p);
    foreach my $g (@p) {
        unless(exists $bg_genes{$g}) { next; }
        unless(exists $ubiq_gene{$g}) {
            $tissue_genes{$tissue}{$g}++;
            $tspfc_genes{$g}++;
        }
    }
}

# Indexing disease-tissue associations
my $desc; my %disease_desc = (); my %disease_tissue = ();
foreach (@dt) {
    @p = split '\t', $_;
    $p[1] =~ s/ /_/g; $p[1] =~ s/\,//g;
    $disease_desc{$disease} = $p[1];
    $disease_tissue{$p[0]} = $p[4];
}

# Indexing disease class and genes
my ($disease, $class, $ngenes, $nt, $nu); my @q;
my %disease_class = ();
my %disease_all_genes = (); my %disease_tspfc_genes = (); my %disease_ubiq_genes = ();
my %class_all_genes = (); my %class_tspfc_genes = (); my %class_ubiq_genes = ();
my %all_genes = (); my %all_tspfc_genes = (); my %all_ubiq_genes = ();
foreach (@dis) {
    @p = split '\t', $_;
    $disease = shift(@p);
    $desc = $disease_desc{$disease};
    #unless(exists $disease_tissue{$disease}) { next; }
    #$tissue = $disease_tissue{$disease};

    @q = split ' ', $p[0];
    $ngenes = pop(@q);
    $class = pop(@q);
    $disease_class{$disease} = $class;
    # $desc = join '_', @q; $desc =~ s/\,//g;
    # $disease_desc{$disease} = $desc;
    shift(@p);

    $nt = 0; $nu = 0;
    foreach my $g (@p) {
        if(exists $tspfc_genes{$g}) {
            $nt++; #print " T:$g";
            $disease_tspfc_genes{$disease}{$g}++;
            $disease_all_genes{$disease}{$g}++;
            unless(($class eq 'Multiple') or ($class eq 'Unclassified')) {
                $class_tspfc_genes{$class}{$g}++; $class_all_genes{$class}{$g}++;
                $all_tspfc_genes{$g}++; $all_genes{$g}++;
            }
        }
        elsif(exists $ubiq_gene{$g}) {
            $nu++; #print " U:$g";
            $disease_ubiq_genes{$disease}{$g}++;
            $disease_all_genes{$disease}{$g}++;
            unless(($class eq 'Multiple') or ($class eq 'Unclassified')) {
                $class_ubiq_genes{$class}{$g}++; $class_all_genes{$class}{$g}++;
                $all_ubiq_genes{$g}++; $all_genes{$g}++;
            }
        }
    }

    if(($nt >= 5) and ($nu >=5)) {
        # if(($nt+$nu) >= 10) {
        print "$disease\t$class\t$ngenes\t$nt\t$nu\t$desc\n"; }
}
#exit;

print "\n";
# For each disease, creating standards, traingin and testing
my %all_pos_genes = (); my %tspfc_pos_genes = (); my %ubiq_pos_genes = (); # +ve genes hash
my ($npos, $nneg); my $min_num_genes = 5; my @rand_genes; # number of genes
# my ($all_pos_gfile, $tspfc_pos_gfile, $ubiq_pos_gfile, $pos_hashref); # files for +ve genes
my $pos_hashref;
my ($all_gfile, $tspfc_gfile, $ubiq_gfile); my $tissue_net; # files for all genes
my ($all_mat, $tspfc_mat, $ubiq_mat);
my ($all_dat, $tspfc_dat, $ubiq_dat);
my ($all_cost, $tspfc_cost, $ubiq_cost);
my ($all_acc, $tspfc_acc, $ubiq_acc);

# my ($aa_mat, $at_mat, $ta_mat, $tt_mat, $ua_mat, $ut_mat, $tempfile); # matrix and temp file
# my ($aa_dat, $at_dat, $ta_dat, $tt_dat, $ua_dat, $ut_dat); # svm data
# my ($aa_cost, $aa_acc, $at_cost, $at_acc, $ta_cost, $ta_acc, $tt_cost, $tt_acc); # svm training
# my ($ua_cost, $ua_acc, $ut_cost, $ut_acc); # svm testing

print TB "#Disease\tDesc\tClass\tTissue\t";
print TB "No.Genes.All\tNo.Genes.Tspfc\tNo.Genes.Ubiq\t";
print TB "Acc.All\tAcc.Tspfc\tAcc.Ubiq\n";

# print TB "AA.acc\tAT.acc\tTA.acc\tTT.acc\tUA.acc\tUT.acc\n";

foreach my $disease (keys %disease_desc) {
    unless(exists $disease_tissue{$disease}) { next; }
    if(scalar(keys %{$disease_tspfc_genes{$disease}}) < 5) { next; }
    $desc = $disease_desc{$disease};
    $class = $disease_class{$disease};
    $tissue_net = $tissuenet_dir.$disease_tissue{$disease}.'.dab';

    # print TB "$disease\t$desc\t$class\t$disease_tissue{$disease}";

    $time = runtime();
    print "$time: $disease\t$desc\n\t$time: printing gene-lists ...\n";

    # $all_pos_gfile = $desc.'_all_pos_genes.txt';
    # $tspfc_pos_gfile = $desc.'_tspfc_pos_genes.txt';
    # $ubiq_pos_gfile = $desc.'_ubiq_pos_genes.txt';

    # open my $APG, ">$all_pos_gfile";
    # open my $TPG, ">$tspfc_pos_gfile";
    # open my $UPG, ">$ubiq_pos_gfile";

    $all_gfile = $desc.'_all_genes.txt';
    $tspfc_gfile = $desc.'_tspfc_genes.txt';
    $ubiq_gfile = $desc.'_ubiq_genes.txt';
    
    open my $AG, ">$all_gfile";
    open my $TG, ">$tspfc_gfile";
    open my $UG, ">$ubiq_gfile";

    # Printing all disease gene lists
    $npos = scalar(keys %{$disease_all_genes{$disease}});
    # $pos_hashref = print_pos($disease_all_genes{$disease}, $APG, $AG);
    $pos_hashref = print_pos($disease_all_genes{$disease}, $AG);
    %all_pos_genes = %{$pos_hashref};
    # print TB "\t$npos";

    @rand_genes = shuffle(keys %all_genes);
    print_neg($npos, \@rand_genes, \%all_pos_genes, \%class_all_genes, $AG);
    close $AG; # close $APG;

    # Printing tspfc disease gene lists
    $npos = scalar(keys %{$disease_tspfc_genes{$disease}});
    # $pos_hashref = print_pos($disease_tspfc_genes{$disease}, $TPG, $TG);
    $pos_hashref = print_pos($disease_tspfc_genes{$disease}, $TG);
    %tspfc_pos_genes = %{$pos_hashref};
    # print TB "\t$npos";

    @rand_genes = shuffle(keys %all_tspfc_genes);
    print_neg($npos, \@rand_genes, \%tspfc_pos_genes, \%class_tspfc_genes, $TG);
    close $TG; # close $TPG;

    # Printing ubiq disease gene lists
    $npos = scalar(keys %{$disease_ubiq_genes{$disease}});
    # $pos_hashref = print_pos($disease_ubiq_genes{$disease}, $UPG, $UG);
    $pos_hashref = print_pos($disease_ubiq_genes{$disease}, $UG);
    %ubiq_pos_genes = %{$pos_hashref};
    # print TB "\t$npos";

    @rand_genes = shuffle(keys %all_ubiq_genes);
    print_neg($npos, \@rand_genes, \%ubiq_pos_genes, \%class_ubiq_genes, $UG);
    close $UG; # close $UPG;

    $time = runtime(); print "\t$time: getting association matrices ...\n";

    # Getting association data from tissue-network and converting them to svm # data
    $all_dat = write_svmdat($tissue_net, $all_gfile, $in_bglist, 'all', \%all_pos_genes);
    $tspfc_dat = write_svmdat($tissue_net, $tspfc_gfile, $in_bglist, 'tspfc', \%tspfc_pos_genes);
    $ubiq_dat = write_svmdat($tissue_net, $ubiq_gfile, $in_bglist, 'ubiq', \%ubiq_pos_genes);
    # $aa_dat = write_svmdat($tissue_net, $all_gfile, $all_pos_gfile, 'aa', \%all_pos_genes);
    # $at_dat = write_svmdat($tissue_net, $all_gfile, $tspfc_pos_gfile, 'at', \%all_pos_genes);
    # $ta_dat = write_svmdat($tissue_net, $tspfc_gfile, $all_pos_gfile, 'ta', \%tspfc_pos_genes);
    # $tt_dat = write_svmdat($tissue_net, $tspfc_gfile, $tspfc_pos_gfile, 'tt', \%tspfc_pos_genes);
    # $ua_dat = write_svmdat($tissue_net, $ubiq_gfile, $all_pos_gfile, 'ua', \%ubiq_pos_genes);
    # $ut_dat = write_svmdat($tissue_net, $ubiq_gfile, $tspfc_pos_gfile, 'ut', \%ubiq_pos_genes);

    $time = runtime(); print "\t$time: svm parameter optimization & prediction ...\n";

    # Linear-SVM grid search
    ($all_cost, $all_acc) = svm_gridsearch($all_dat);
    ($tspfc_cost, $tspfc_acc) = svm_gridsearch($tspfc_dat);
    ($ubiq_cost, $ubiq_acc) = svm_gridsearch($ubiq_dat);
    # ($aa_cost, $aa_acc) = svm_gridsearch($aa_dat);
    # ($at_cost, $at_acc) = svm_gridsearch($at_dat);
    # ($ta_cost, $ta_acc) = svm_gridsearch($ta_dat);
    # ($tt_cost, $tt_acc) = svm_gridsearch($tt_dat);

    # Linear-SVM training and testing
    # $ua_acc = svm_traintest($ta_dat, $ta_cost, $ua_dat);
    # $ut_acc = svm_traintest($tt_dat, $tt_cost, $ut_dat);

    print TB "\t$all_acc\t$tspfc_acc\t$ubiq_acc\n";
    # print TB "\t$aa_acc\t$at_acc\t$ta_acc\t$tt_acc\t$ua_acc\t$ut_acc\n";
}

close TB;
$time = runtime(); print "\n$time: DONE\n\n";

