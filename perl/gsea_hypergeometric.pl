#!/usr/bin/perl
use strict;
use warnings;
use lib '/Genomics/Users/arjunk/software/lib/perl5/';
use PDF;
use Getopt::Long;
use Pod::Usage;
use Time::SoFar qw(runtime);
use Graph::Undirected;

sub index_list;
sub gene2fam;
sub parse_gmt;
sub bring2_commonbg;
sub updatems;
sub get_esfdr;


my ($help, $igmt1, $igmt2, $ibg1, $ibg2, $igsfilt, $orich, $ogene, $ouniq, $ibest, $time, @ref, @p, $itrim, $iortho);
my $opar = 'tab'; my $omats = 'lor'; my $omatc;
my $iming = 5; my $imaxg1 = my $imaxg2 = 200;
my $imincom = 3; my $inoident = 1;
# my $imincom = 0; my $inoident = 0;
my $iqval = 0.01; my ($iovlp, $ijac, $ilor, $ipval);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions( 'help' => \$help,
         'igmt1=s' => \$igmt1,
         'igmt2=s' => \$igmt2,
        'iortho=s' => \$iortho,
          'ibg1=s' => \$ibg1,
          'ibg2=s' => \$ibg2,
       'igsfilt=s' => \$igsfilt,
         'iming=i' => \$iming,
        'imaxg1=i' => \$imaxg1,
        'imaxg2=i' => \$imaxg2,
       'imincom=i' => \$imincom,
      'inoident=i' => \$inoident,
         'iqval=f' => \$iqval,
         'iovlp=f' => \$iovlp,
          'ijac=f' => \$ijac,
          'ilor=f' => \$ilor,
         'ipval=f' => \$ipval,
         'ibest=s' => \$ibest,
          'opar=s' => \$opar,
         'omats=s' => \$omats,
           'omatc' => \$omatc,
           'itrim' => \$itrim,
         'orich=s' => \$orich,
         'ogene=s' => \$ogene,
         'ouniq=s' => \$ouniq) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if($iming < $imincom) { $iming = $imincom; }

unless($orich) {
    (my $tag1 = $igmt1) =~ s/.*\///g; $tag1 =~ s/\.[a-z]*$//g;
    (my $tag2 = $igmt2) =~ s/.*\///g; $tag2 =~ s/\.[a-z]*$//g;
    $orich = $tag1.'-'.$tag2.'.gsea.out';
    print "\nyour output is stored in $orich\n"; }

my %gsym = my %gname = ();
if($ogene) {
    my $isym = '/Genomics/ogtr04/arjunk/data/mappings/human_gene-info_ncbi.txt';
    open GH, "$isym" or die "Can't open $isym!";
    while (<GH>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $gsym{$p[0]} = $p[1];
        $gname{$p[0]} = $p[2]; }
    close GH; }

if($ibest) {
    if($opar eq 'mat') {
        print "\nwarning: for recording only best-matches, a table is preferred.\n";
        $opar = 'tab'; }
    undef $iovlp; undef $ijac; undef $ilor; undef $ipval; }


$time = runtime(); print "\n$time: Assigning genes-to-genesets and trimming them to bg ...";
my %gene_fam = ();
if($iortho) {
    %gene_fam = gene2fam($iortho);
    @ref = parse_gmt($igmt1, $imaxg1, \%gene_fam); }
else {
    @ref = parse_gmt($igmt1, $imaxg1); }

my %gs1_genes = %{$ref[0]};
my %gs1_desc = %{$ref[1]};
my %gs1_allg = %{$ref[2]};

unless($igmt2) {
    $igmt2 = $igmt1;
    $imaxg2 = $imaxg1; }
if($iortho) {
    @ref = parse_gmt($igmt2, $imaxg2, \%gene_fam); }
else {
    @ref = parse_gmt($igmt2, $imaxg2); }

my %gs2_genes = %{$ref[0]};
my %gs2_desc = %{$ref[1]};
my %gs2_allg = %{$ref[2]};

my %bg1 = my %bg2 = my %all_genes = ();
my %gs1_size = my %gs2_size = ();

if($ibg2) {
    %bg2 = index_list($ibg2);
    my($size, $allg) = bring2_commonbg(\%gs1_genes, \%bg2, $imaxg1);
    %gs1_size = %{$size};
    map { $all_genes{$_}++; } keys %{$allg}; }
else {
    my($size, $allg) = bring2_commonbg(\%gs1_genes, \%gs2_allg, $imaxg1);
    %gs1_size = %{$size};
    map { $all_genes{$_}++; } keys %{$allg}; }

if($ibg1) {
    %bg1 = index_list($ibg1);
    my($size, $allg) = bring2_commonbg(\%gs2_genes, \%bg1, $imaxg2);
    %gs2_size = %{$size};
    map { $all_genes{$_}++; } keys %{$allg}; }
else {
    my($size, $allg) = bring2_commonbg(\%gs2_genes, \%gs1_allg, $imaxg2);
    %gs2_size = %{$size};
    map { $all_genes{$_}++; } keys %{$allg}; }


my $univ = scalar keys %all_genes;

my @gs1_array = sort keys %gs1_genes;
my @gs2_array = sort keys %gs2_genes;

if($igsfilt) {
    my %filt = index_list($igsfilt);

    @gs1_array = grep { exists $filt{$_} } @gs1_array;
    @gs2_array = grep { exists $filt{$_} } @gs2_array; }

print "\n\tGS1: $igmt1\n\tNo. genesets: ", scalar @gs1_array, "\n\t";
print "No. genes: ", scalar keys %gs1_allg, "\n\n";
print "\tGS2: $igmt2\n\tNo. genesets: ", scalar @gs2_array, "\n\t";
print "No. genes: ", scalar keys %gs2_allg, "\n\n";
print "\tBackground: ", scalar keys %all_genes, "\n\n";


my %gs1_count = my %gs1_mean = my %gs1_sd = ();
my %gs2_count = my %gs2_mean = my %gs2_sd = ();
my %gspair_score = ();

if($omatc) {
    foreach (@gs1_array) {
        $gs1_count{$_} = $gs1_mean{$_} = $gs1_sd{$_} = 0; } }


open OUT, ">$orich";

if($opar eq 'mat') {
    print OUT "\t";
    foreach (@gs2_array) {
        $gs2_count{$_} = $gs2_mean{$_} = $gs2_sd{$_} = 0;
        print OUT "\t$_"; } print OUT "\n\t";
    foreach (@gs2_array) { print OUT "\t$gs2_desc{$_} ($gs2_size{$_})"; }
    print OUT "\n"; }
else {
    print OUT "#gs1\tgs1.ngenes\tgs1.desc\t";
    print OUT "gs2\tgs2.ngenes\tgs2.desc\t";
    print OUT "ncommon\tovlp\tjac\tlor\tpvalue\tes.qvalue\n"; }

$time = runtime(); print "$time: Starting pair-wise comparisons ...";
my $count = 0;
my %num_common  = my %com_genes = ();
my %pvalue = my %lor = my %overlap = my %jaccard = my %allp = ();
my ($best_ncom, $best_ovlp, $best_jac, $best_pval, $best_lor, $best_match);

GS1: foreach my $gs1 (@gs1_array) {
    my $size1 = $gs1_size{$gs1};
    if($size1 == 0) { next GS1; }
    $count++;

    if($ibest) {
        $best_match = '';
        $best_ncom = 0;
        $best_ovlp = 0;
        $best_jac = 0;
        $best_pval = 1;
        $best_lor = 0; }

    if($opar eq 'mat') {
        unless($omatc) {
            print OUT "$gs1\t$gs1_desc{$gs1} ($gs1_size{$gs1})"; } }

    GS2: foreach my $gs2 (@gs2_array) {
        if($inoident and ($gs1 eq $gs2)) {
            if($opar eq 'mat') { unless($omatc) { print OUT "\tNA"; } }
            next GS2; }

        my $tag;
        if($igmt2 eq $igmt1) { $tag = join '__', sort($gs1, $gs2); }
        else { $tag = join '__', ($gs1, $gs2); }

        my $size2 = $gs2_size{$gs2}; if($size2 == 0) { next GS2; }
        my $min = $size2; if($size1 < $size2) { $min = $size1; }

        my @comg = grep { exists $gs1_genes{$gs1}{$_} } keys %{$gs2_genes{$gs2}};
        my $ncom = scalar @comg;
        #if(($opar ne 'mat') and ($ncom < $imincom)) {
        #    next GS2; }
        foreach my $g (@comg) { $com_genes{$tag}{$g}++; }

        my %tempu = %{$gs1_genes{$gs1}};
        @tempu{keys %{$gs2_genes{$gs2}}} = values %{$gs2_genes{$gs2}};
        my $nuni = scalar keys %tempu;

        my $ovlp = sprintf("%.6f", ($ncom/$min));
        my $jac = sprintf("%.6f", ($ncom/$nuni));
        my $pval = hypergeometric_tail($univ, $size1, $size2, $ncom);
        my $lor = 0; if($ncom > 0) {
            $lor = sprintf("%.6f", log(($ncom/$size1)/($size2/$univ))/log(2)); }

        if($opar eq 'mat') {
            my $mats = 0;
            if($omats eq 'ncom') {
                $mats = $ncom; }
            elsif($omats eq 'ovlp') {
                $mats = $ovlp; }
            elsif($omats eq 'jac') {
                $mats = $jac; }
            elsif($omats eq 'pval') {
                if($pval==0) { $mats = 325 }
                else {
                    $mats = -1*log($pval)/log(10); } }
            elsif($omats eq 'lor') {
                $mats = $lor; }

            if($omatc) {
                $gspair_score{$gs1}{$gs2} = $mats;

                $gs1_count{$gs1}++;
                my $om = $gs1_mean{$gs1};
                $gs1_mean{$gs1} += (($mats - $om) / $gs1_count{$gs1});
                $gs1_sd{$gs1} += ($mats - $om)*($mats - $gs1_mean{$gs1});

                $gs2_count{$gs2}++;
                $om = $gs2_mean{$gs2};
                $gs2_mean{$gs2} += (($mats - $om) / $gs2_count{$gs2});
                $gs2_sd{$gs2} += ($mats - $om)*($mats - $gs2_mean{$gs2}); }
            else {
                print OUT "\t$mats";
                #if($omats eq 'lor') {
                #    if($mats > 0) {
                #        if($ncom < $imincom) { print OUT "\tNA"; } }
                #    else { print OUT "\t$mats"; } }
                #else {
                #    if($ncom < $imincom) { print OUT "\tNA"; }
                #    else { print OUT "\t$mats"; } }
            }
            next GS2; }

        if($ibest) {
            if($ncom < $imincom) { next GS2; }
            if((($ibest eq 'ovlp') and ($ovlp > $best_ovlp))
                    or (($ibest eq 'jac') and ($jac > $best_jac))
                    or (($ibest eq 'pval') and ($pval < $best_pval))
                    or (($ibest eq 'lor') and ($lor > $best_lor))) {
                $best_ncom = $ncom;
                $best_ovlp = $ovlp;
                $best_jac = $jac;
                $best_pval = $pval;
                $best_lor = $lor;
                $best_match = $gs2; }
            else { next GS2; } }
        else {
            if($iovlp) { $allp{$tag} = $ovlp;
                if(($ncom < $imincom) or ($ovlp < $iovlp)) { next GS2; } }
            elsif($ijac) { $allp{$tag} = $jac;
                if(($ncom < $imincom) or ($jac < $ijac)) { next GS2; } }
            elsif($ilor) { $allp{$tag} = $lor;
                if(($ncom < $imincom) or ($lor < $ilor)) { next GS2; } }
            elsif($ipval) { $allp{$tag} = $pval;
                if(($ncom < $imincom) or ($pval > $ipval)) { next GS2; } }
            elsif($iqval) {
                if($ncom < $imincom) { next GS2; } }

            $num_common{$tag} = $ncom;
            $overlap{$tag} = $ovlp;
            $jaccard{$tag} = $jac;
            $pvalue{$tag} = $pval;
            $lor{$tag} = $lor; } }

    if($opar eq 'mat') { unless($omatc) { print OUT "\n"; next GS1; } }

    if($ibest) {
        unless($best_match eq '') {
            my $tag = $gs1.'__'.$best_match;
            $num_common{$tag} = $best_ncom;
            $overlap{$tag} = $best_ovlp;
            $jaccard{$tag} = $best_jac;
            $lor{$tag} = $best_lor;
            $pvalue{$tag} = $best_pval; } }

    if(($count % 500) == 0) {
        $time = runtime(); print "\n\t$time: $count ..."; } }

if($opar eq 'mat') {
    if($omatc) {
        foreach (keys %gs1_sd) {
            $gs1_sd{$_} = sqrt($gs1_sd{$_} / ($gs1_count{$_}-1));
            if($gs1_sd{$_} == 0) { die "\nsd is zero for $_ !\n\n"; } }

        foreach (keys %gs2_sd) {
            $gs2_sd{$_} = sqrt($gs2_sd{$_} / ($gs2_count{$_}-1));
            if($gs2_sd{$_} == 0) { die "\nsd is zero for $_ !\n\n"; } }

        (my $otab = $orich) =~ s/\.mat$/\.txt/g;
        open TAB, ">$otab";
        print TAB "#gs1\tgs1.desc\tgs1.ngenes\tgs2\tgs2.desc\tgs2.ngenes\toverlap\tzscore\n";

        $time = runtime(); print "\n$time: Background-correcting overlap measures ...";

        foreach my $gs1 (@gs1_array) {
            print OUT "$gs1\t$gs1_desc{$gs1} ($gs1_size{$gs1})";

            foreach my $gs2 (@gs2_array) {
                if($inoident and ($gs1 eq $gs2)) {
                    print OUT "\tNA"; next; }

                my $s = $gspair_score{$gs1}{$gs2};
                my $z1 = ($s - $gs1_mean{$gs1}) / $gs1_sd{$gs1};
                my $z2 = ($s - $gs2_mean{$gs2}) / $gs2_sd{$gs2};
                my $z = ($z1 + $z2) / sqrt(2);
                print OUT "\t$z";

                print TAB "$gs1\t$gs1_desc{$gs1}\t$gs1_size{$gs1}\t";
                print TAB "$gs2\t$gs2_desc{$gs2}\t$gs2_size{$gs2}\t";
                print TAB "\t$gspair_score{$gs1}{$gs2}\t$z\n";
            }
            print OUT "\n"; }
        close TAB; }
    close OUT;

#     if($itrim) {
#         $time = runtime(); print "\n$time: Trimming LOR matrix ...";
#         my $mat = $orich;
#         (my $tab = $orich) =~ s/\.mat/\.txt/g;
#         (my $out = $orich) =~ s/\.mat/_trimmed\.mat/g;
#         `trim_lor-matrix.pl $mat $tab $out`; }

    $time = runtime(); print "\n$time: DONE\n\n";
    exit; }

$time = runtime(); print "\n\n$time: Performing multiple-hypothesis correction ...";
my %esq = get_esfdr(\%pvalue);

$time = runtime(); print "\n\n$time: Printing results ...";

my %print_set;
if($iovlp or ($ibest and ($ibest eq 'ovlp'))) { %print_set = %overlap; }
elsif($ijac or ($ibest and ($ibest eq 'jac'))) { %print_set = %jaccard; }
elsif($ilor or ($ibest and ($ibest eq 'lor'))) { %print_set = %lor; }
elsif($ipval or ($ibest and ($ibest eq 'pval'))) { %print_set = %pvalue; }
elsif($iqval or ($ibest and ($ibest eq 'qval'))) { %print_set = %esq; $iqval = (-1)*log($iqval)/log(10); }

my ($gs1, $gs2);
foreach my $tag (sort {$print_set{$b} <=> $print_set{$a}} keys %print_set) {
    unless($iovlp or $ijac or $ilor or $ipval or $ibest) {
        if($print_set{$tag} < $iqval) { last; } }

    ($gs1, $gs2) = split '__', $tag;

    if($inoident) {
        print OUT "$gs1\t$gs1_size{$gs1}\t$gs1_desc{$gs1}\t";
        print OUT "$gs2\t$gs2_size{$gs2}\t$gs2_desc{$gs2}\t";
        print OUT "$num_common{$tag}\t$overlap{$tag}\t";
        print OUT "$jaccard{$tag}\t$lor{$tag}\t";
        print OUT sprintf("%.6g\t%.6g\n", $pvalue{$tag}, $esq{$tag}); }
    else {
        unless($gs1 eq $gs2) { next; }
        print OUT "$gs1\t$gs1_desc{$gs1}\t$gs1_size{$gs1}\t$gs2_size{$gs2}\t";
        print OUT "$num_common{$tag}\t$overlap{$tag}\t";
        print OUT "$jaccard{$tag}\t$lor{$tag}\t";
        print OUT sprintf("%.6g\t%.6g\n", $pvalue{$tag}, $esq{$tag}); } }
close OUT;

if(($opar eq 'tab') and $ogene) {
    open OUTG, ">$ogene";
    if($inoident) {
        print OUTG "#gs1\tgs1.desc\tgs2\tgs2.desc\tgene\tsymbol\tdesc\n";
        foreach my $tag (sort {$print_set{$b} <=> $print_set{$a}} keys %print_set) {
            unless($iovlp or $ijac or $ilor or $ipval or $ibest) {
                if($print_set{$tag} < $iqval) { last; } }

            my ($gs1, $gs2) = split '__', $tag;

            foreach my $g (keys %{$com_genes{$tag}}) {
                my $y = my $d = '--';
                if(exists $gsym{$g}) {
                    $y = $gsym{$g}; $d = $gname{$g}; }
                print OUTG "$gs1\t$gs1_desc{$gs1}\t$gs2\t$gs2_desc{$gs2}\t$g\t$y\t$d\n"; } } }
    else {
        print OUTG "#gs\tgs.desc\tsp1.genes\tsp2.genes\tcom.genes\n";
        foreach my $tag (sort {$print_set{$b} <=> $print_set{$a}} keys %print_set) {
            unless($iovlp or $ijac or $ilor or $ipval or $ibest) {
                if($print_set{$tag} < $iqval) { last; } }

            my ($gs1, $gs2) = split '__', $tag;
            unless($gs1 eq $gs2) { next; }

            print OUTG "$gs1\t$gs1_desc{$gs1}";
            print OUTG "\t", scalar grep { not exists $com_genes{$tag}{$_} } keys %{$gs1_genes{$gs1}};
            print OUTG "\t", scalar grep { not exists $com_genes{$tag}{$_} } keys %{$gs2_genes{$gs1}};
            print OUTG "\t", scalar keys %{$com_genes{$tag}};
            my $str = join ';', grep { not exists $com_genes{$tag}{$_} } keys %{$gs1_genes{$gs1}}; if(length($str)==0) { $str = '--'; }
            print OUTG "\t", $str;
            $str = join ';', grep { not exists $com_genes{$tag}{$_} } keys %{$gs2_genes{$gs1}}; if(length($str)==0) { $str = '--'; }
            print OUTG "\t", $str;
            $str = join ';', keys %{$com_genes{$tag}}; if(length($str)==0) { $str = '--'; }
            print OUTG "\t", $str, "\n"; } }
    close OUTG; }

if(($opar eq 'tab') and $ouniq) {
    open OUTU, ">$ouniq";
    # print OUTU "#GS\tGS.Desc\tGS.Size\tGS.Ovlp\n";

    my %all_gs_desc = my %all_gs_size = ();
    foreach my $gs (@gs1_array) {
        $all_gs_desc{$gs} = $gs1_desc{$gs};
        $all_gs_size{$gs} = $gs1_size{$gs}; }
    if($igmt2 ne $igmt1) {
        foreach my $gs (@gs2_array) {
            $all_gs_desc{$gs} = $gs2_desc{$gs};
            $all_gs_size{$gs} = $gs2_size{$gs}; } }

    my $graph = Graph::Undirected->new;
    foreach my $tag (sort {$print_set{$b} <=> $print_set{$a}} keys %print_set) {
        unless($iovlp or $ijac or $ilor or $ipval or $ibest) {
            if($print_set{$tag} < $iqval) { last; } }

        my ($gs1, $gs2) = split '__', $tag;
        # $graph->add_weighted_edge($gs1, $gs2, $print_set{$tag}); }
        $graph->add_edge($gs1, $gs2); }

    foreach my $gs (keys %all_gs_desc) {
        unless($graph->has_vertex($gs)) {
            $graph->add_vertex($gs); } }

    my @cc = $graph->connected_components();

    print "\n";
    foreach my $c (@cc) {
        print "\n";
        foreach my $gs (@{$c}) {
            print "$gs|$all_gs_desc{$gs}|$all_gs_size{$gs}\n"; }

        my $ngs = scalar @{$c};
        if($ngs == 1) {
            my $gs = ${$c}[0];
            # print OUTU "$gs\t$all_gs_desc{$gs}\t$all_gs_size{$gs}\tNA\n";
            print OUTU "$gs\t$all_gs_desc{$gs} ($all_gs_size{$gs})\t";
            my @gsg = ();
            if(exists $gs1_genes{$gs}) { @gsg = keys %{$gs1_genes{$gs}}; }
            elsif(exists $gs2_genes{$gs}) { @gsg = keys %{$gs2_genes{$gs}}; }
            print OUTU join "\t", @gsg; print OUTU "\n";
            next; }

        my %rel_pairs = my %gs_ovlp = ();
        for(my $i=0; $i<$#{$c}; $i++) {
            my $gs1 = ${$c}[$i];
            for(my $j=($i+1); $j<=$#{$c}; $j++) {
                my $gs2 = ${$c}[$j];

                my $tag = join '__', sort($gs1, $gs2);
                unless(exists $allp{$tag}) { print "\n$tag : no score\n"; exit; }
                if(exists $allp{$tag}) {
                    $gs_ovlp{$gs1} += $allp{$tag};
                    $gs_ovlp{$gs2} += $allp{$tag}; }
        # $gs_ovlp{$gs1} += $graph->get_edge_weight($gs1, $gs2);
        # $gs_ovlp{$gs2} += $graph->get_edge_weight($gs1, $gs2); } }

                unless($graph->has_edge($gs1, $gs2)) { next; }
                $rel_pairs{$gs1}{$gs2}++; $rel_pairs{$gs2}{$gs1}++; } }

        # my @tempa = sort { $all_gs_size{$b} <=> $all_gs_size{$a} } @{$c};
        # $ugs = $tempa[0];
        # my @tempo = (); map { $tempo{$_}++; } values %gs_ovlp;

        my @tempa = sort { $gs_ovlp{$b} <=> $gs_ovlp{$a} } @{$c};
        my %uniq_gs = my %redu_gs = ();
        foreach my $gs (@tempa) {
            if(exists $redu_gs{$gs}) { next; }
            $uniq_gs{$gs}++;

            # print OUTU "$gs\t$all_gs_desc{$gs}\t$all_gs_size{$gs}\t$gs_ovlp{$gs}\n";
            print OUTU "$gs\t$all_gs_desc{$gs} ($all_gs_size{$gs})\t";
            my @gsg = ();
            if(exists $gs1_genes{$gs}) { @gsg = keys %{$gs1_genes{$gs}}; }
            elsif(exists $gs2_genes{$gs}) { @gsg = keys %{$gs2_genes{$gs}}; }
            print OUTU join "\t", @gsg; print OUTU "\n";

            foreach my $ogs (keys %{$rel_pairs{$gs}}) {
                $redu_gs{$ogs}++; } } }

        # print OUTU "$ugs\t$all_gs_desc{$ugs}\t$all_gs_size{$ugs}\n"; }
    close OUTU; }
 
$time = runtime();
print "\n\n$time: DONE\n\n";


# Index list
sub index_list {
    my $list = shift;
    my %index = ();

    open GH, "$list" or die "Can't open $list!";
    while (<GH>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        $index{$p[0]}++; }
    close GH;

    return %index; }


# Assigns to family to genes
sub gene2fam {
    my $ortho = shift;
    my %gfam = ();

    open FH, "$ortho" or die "Can't open $ortho!";
    while (<FH>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        my $fam = shift @p;
        foreach my $spg (@p) {
            (my $g = $spg) =~ s/^[a-zA-Z]*\|//g;
            $gfam{$g}{$fam}++; } }
    close FH;

    return %gfam; }


# Creates a gs -> genes hash from a .gmt file, or a gene -> gs hash otherwise
sub parse_gmt {
    my $gmt = shift;
    my $maxg = shift;
    my $gfam = shift;

    my %genes = my %desc = my %allg = ();
    
    open GMT, "$gmt" or die "Can't open $gmt!";
    while (<GMT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;

        my $gs = shift @p;
        ($desc{$gs} = $p[0]) =~ s/ \([0-9]*\)$//g;

        shift @p;
        foreach my $g (@p) {
            if($iortho) {
                unless(exists $gfam->{$g}) { next; }
                foreach my $f (keys %{$gfam->{$g}}) {
                    $genes{$gs}{$f}++;
                    $allg{$f}++; } }
            else {
                $genes{$gs}{$g}++;
                $allg{$g}++; } }

        if((scalar keys %{$genes{$gs}} > $maxg) or
            (scalar keys %{$genes{$gs}} < $iming)) {
            delete $genes{$gs}; } }
    close GMT;

    return (\%genes, \%desc, \%allg); }


# Restrict genes in geneset to specified background
sub bring2_commonbg {
    my $gsref = shift;
    my $generef = shift;
    my $maxg = shift;

    my %size = my %allg = ();

    foreach my $gs (keys %{$gsref}) {
        $size{$gs} = 0;
        foreach my $g (keys %{$gsref->{$gs}}) {
            if(exists $generef->{$g}) {
                $allg{$g}++;
                $size{$gs}++; }
            else {
                delete $gsref->{$gs}->{$g}; } } }

    foreach my $gs (keys %size) {
        if(($size{$gs} < $iming) or ($size{$gs} > $maxg)) {
            delete $gsref->{$gs}; } }

    return (\%size, \%allg); }


# Update mean and sd
sub updatems {
    my $val = shift;
    my $mean = shift;
    my $sd = shift;
    my $n = shift;

    my $ume = $mean + (($val - $mean) / $n);
    my $usd = $sd + ($val - $mean)*($val - $ume);

    return ($ume, $usd); }


# Performs multiple-hypothesis correction using BH method
sub get_esfdr {
    my $pref = shift;
    
    my $minp = 1; my $ntests = 0;
    foreach (keys %{$pref}) {
        $ntests++;
        if(($pref->{$_} != 0) 
                and ($pref->{$_} < $minp)) {
            $minp = $pref->{$_}; } }

    my %esfdr = (); my $rank = $ntests;
    foreach (sort {$pref->{$b} <=> $pref->{$a}} keys %{$pref}) {
        if($pref->{$_} == 0) {
            $esfdr{$_} = (-1)*log( ($minp * $ntests) / ($rank * 10)) / log(10); }
        else {
            $esfdr{$_} = (-1)*log(($pref->{$_} * $ntests) / $rank)/log(10); }

        if($esfdr{$_} < 0) { $esfdr{$_} = 0; }

        $rank--; }

    return %esfdr; }

__END__

=head1

Calculates gene-overlap between all pairs of genesets using the hypergeometric
test.

=head1 USAGE

./gsea_hypergeometric.pl [--igmt1 GENESET1_GMT] [--igmt2 GENESET2_GMT]
[--iortho GENE-FAMLIES] [--iming MIN_NUM_GENES] [--imaxg1 MAX_NUM_GENES_IN_GS1]
[--imaxg2 MAX_NUM_GENES_IN_GS2] [--iqval QVALUE_CUTOFF] [--iovlp OVERLAP_CUTOFF]
[--ijac JACCARD_CUTOFF] [--ilor LOGODDSRATIO_CUTOFF] [--opar OUTPUT_PAR]
[--orich OUTPUT_FILE] [--ogene OVERLAPPING_GENES] [--ouniq NON-OVERLAPPING_GS]
[--omatc SCORE_CORRECTION] [--help]

=head1 DESCRIPTION

This script takes in two collections of genesets and calculates gene-overlap
between all pairs of genesets across the collections and report enrichment
statistics. The output can be a plain table listing all pairs of genesets, or a
matrix with genesets along the rows & columns, with each cell recording the
log-odds-ratio of the overlap.

=head1 ARGUMENTS

=over 12

=item C<--igmt1>

First collection of genestes in GMT format.

=item C<--igmt2>

(Optional) Second collection of genestes in GMT format. If provided, genesets
across the collections will be compared. If not, genesets within the first
collection will be compared.

=item C<--iortho>

(Optional) This file containing gene-famnily-to-gene mapping can be used to
compare geneset collections from two different species. Each line in this file
should begin with a family_id followed by its members in a tab-delimited format.

=item C<--ibg1>

(Optional) Background genes for collection 1.

=item C<--ibg2>

(Optional) Background genes for collection 2.

=item C<--igsfilt>

(Optional) List of genesets to consider. Backgrounds are built based on the
entire collections, btu comparisons are performed only between genesets that are
present in this list.

=item C<--iming>

(Optional) Minimum no. of genes in a geneset. Genesets with #genes < iming will
be excluded from the analysis. Default is 5.

=item C<--imaxg1>

(Optional) Maximum no. of genes in a geneset in collection 1. Genesets with
#genes > imaxg1 will be excluded from the analysis. Default is 200.

=item C<--imaxg2>

(Optional) Maximum no. of genes in a geneset in collection 2. Genesets with
#genes > imaxg2 will be excluded from the analysis. Default is 200.

=item C<--imincom>

(Optional) Minimum no. of genes that should be common between a pair of genesets
to be considered further for analysis or printing.

=item C<--inoident>

(Optional) 1: genesets with identical IDs will not be considered for further
analysis/printing. 0: identical genesets will also be taken seriously. Default
is 1.

=item C<--iqval>

(Optional) Q-value (FDR) cut-off to declare a geneset pair significant. Default
is 0.01. --iovlp, --ijac or --ilor, if provided, take precedence.

=item C<--iovlp>

(Optional) Overlap cut-off. Overlap = |A int B| / min(|A|, |B|). Ranges between
0 and 1. --iqval is ignored.

=item C<--ijac>

(Optional) Jaccard cut-off. Jaccard = |A int B| / |A uni B|. Ranges between
0 and 1. --iqval is ignored.

=item C<--ilor>

(Optional) Log-Odds-Ratio cut-off. LOR = log_2[(|A int B|/|A|)/(|B|/Total)].
Over-presentation among sets will yield +ve values and under-representations
will yield -ve values, both unbounded. --iqval is ignored.

=item C<--ipval>

(Optional) P-value cut-off to declare a geneset pair significant. --iqval is
ignored.

=item C<--ibest>

(Optional) Only best-matches from the second collection to each geneset in the
first will be reported in a table based on the argument provided - 'ovlp',
'jac', 'pval' or 'lor'. --iqval, --iovlp or --ijac, if set, are ignored. All
matches with >=3 genes in common are considered.

=item C<--opar>

(Optional) Output file format/mode. Default is 'tab'.

'tab' : all pairs of genesets with overlap Q-value < iqval should be listed
along with their overlap statistics.

'mat' : the genesets from the two collections should be placed along the rows &
columns of a matrix, with each cell recording the measure of the overlap, which
can be chosen using the --omats option.

=item C<--omats>

(Optional) The overlap measure to be printed in the matrix. Can be one of
'ncom', 'ovlp', 'jac' or 'lor'. Should be specified along with --opar set to
'mat'. Default is 'lor'.

=item C<--omatc>

(Optional) If specified along with --opar set to 'mat', the matrix of overlap
measures are background corrected based on the distribution of the scores for
the involved genesets with all other sets.

=item C<--orich>

(Optional) Output file recording output as specified by --opar.

=item C<--ogene>

(Optional) Output file containing genes per siginificant overlap.

=item C<--ouniq>

(Optional) Output file containing non-overlapping genesets. This option makes a lot of
sense only when comparing genesets within a collection.

=item C<--help>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut


