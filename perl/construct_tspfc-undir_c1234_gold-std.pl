#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Graph;
use List::Util qw(shuffle);
use Time::SoFar qw(runtime);

# Subroutines
sub parse_obo;

# Directories
my $exp_dir = '/home/arjunk/data/gene-expression/';

# Input files
# my $itis_gmt = $exp_dir.'hprd/hprd-brenda/human_tissue-specific_gene-expression_pruned-direct-10.gmt';
my $iunp_tis_gmt = $exp_dir.'hprd/hprd-brenda/human_tissue-specific_gene-expression_unprop.gmt';
my $iubq_glist = $exp_dir.'ubiquitous_genes/human/alldb_ubiquitous-genes_human.txt';
my $ibg_glist = '/Genomics/ogtr04/arjunk/projects/human-tissue-fln/standards/human_standard.genes';
my $iobo = '/home/arjunk/data/ontologies/brenda/BrendaTissue.obo';

my ($help, $ipos, $ineg, $ispfc, $inouue, $ifilt); my $ionto = 'genes';
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV < 2);
GetOptions( 'help' => \$help,
          'ipos=s' => \$ipos,
          'ineg=s' => \$ineg,
  'iunp_tis_gmt=s' => \$iunp_tis_gmt,
         'ifilt=s' => \$ifilt,
    'iubq_glist=s' => \$iubq_glist,
     'ibg_glist=s' => \$ibg_glist,
          'iobo=s' => \$iobo,
         'ionto=s' => \$ionto,
          'inouue' => \$inouue,
           'ispfc' => \$ispfc   ) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

if($inouue) { undef $ispfc; }
my ($time, @p, @r, $out_gsdat, $out_gsdab, $e);


# Output files
open STAT, ">human-tissue-standards.txt";


$time = runtime(); print "\n$time: Parsing background and ubiquitous gene lists ...\n";

# Indexing background genes
open BG, "$ibg_glist" or die "Can't open $ibg_glist!";
my %bg_genes = ();
while (<BG>) {
    chomp($_); $bg_genes{$_}++; }
close BG;


my %ubiq_genes = (); my %all_genes = ();
# Indexing ubiquitously expresses genes
open UB, "$iubq_glist" or die "Can't open $iubq_glist!";
while (<UB>) {
    chomp($_);
    unless(exists $bg_genes{$_}) { next; }
    $ubiq_genes{$_}++;
    unless($inouue) { $all_genes{$_}++; } }
close UB;


$time = runtime(); print "\n$time: Parsing tissue-gene annotations ...\n";

# Indexing unpropagated non-ubiquitous tissue-specific genes
open UTS, "$iunp_tis_gmt" or die "Can't open $iunp_tis_gmt!";
my ($t, $desc); my %all_tspfc_genes = (); my %tis_genes = ();

while (<UTS>) {
    chomp($_); @p = split '\t', $_;
    $t = shift(@p); shift(@p);
    %{$tis_genes{$t}} = ();
    foreach my $g (@p) {
        unless(exists $bg_genes{$g}) { next; }
        unless($inouue) { if(exists $ubiq_genes{$g}) { next; } }
        $tis_genes{$t}{$g}++;
        $all_tspfc_genes{$g}++;
        $all_genes{$g}++; } }
close UTS;

my $num_all_genes = scalar keys %all_genes;

my ($num_ubiq_genes, $num_tspfc_genes);
if($inouue) {
    $num_ubiq_genes = $num_tspfc_genes = 0;
    foreach my $g (keys %all_genes) {
        if(exists $ubiq_genes{$g}) { $num_ubiq_genes++; }
        else { $num_tspfc_genes++; } } }
else {
    $num_ubiq_genes = scalar keys %ubiq_genes;
    $num_tspfc_genes = scalar keys %all_tspfc_genes; }

print "\tNo. all genes: $num_all_genes\n";
print "\tNo. ubiq genes: $num_ubiq_genes\n";
print "\tNo. tspfc genes: $num_tspfc_genes\n";


# Parsing ontology
my $root = 'BTO_0000000';
my (@d, @a); my %tis_des = my %tis_anc = ();

$time = runtime(); print "\n$time: Parsing tissue ontology ...\n";
my $obog = parse_obo($iobo);

foreach my $t (keys %tis_genes) {
    @d = $obog->all_predecessors($t);
    @d = grep { exists $tis_genes{$_} } @d; if(scalar @d == 0) { next; }
    foreach (@d) { unless($_ eq $root) { $tis_des{$t}{$_}++; } }

    @a = $obog->all_successors($t);
    @a = grep { exists $tis_genes{$_} } @a; if(scalar @a == 0) { next; }
    foreach (@a) { unless($_ eq $root) { $tis_anc{$t}{$_}++; } } }


# Indexing propagated tissue-specific genes
my $tis_exp_dat = 'human_tissue-exp.dat';
my %tis_prop_genes = (); my %tis_ign_genes = ();
if($ionto eq 'genes') {
    $time = runtime(); print "\n$time: Propagating tissue genes ...\n";

    unless(-e $tis_exp_dat) {
        open TEXP, ">tmp_tis-exp.dat"; }

    foreach my $t (keys %tis_genes) {
        foreach my $g (keys %{$tis_genes{$t}}) {
            $tis_prop_genes{$t}{$g}++;
            $tis_ign_genes{$t}{$g}++; }

        foreach my $dt (keys %{$tis_des{$t}}) {
            unless(exists $tis_genes{$dt}) { next; }
            foreach my $g (keys %{$tis_genes{$dt}}) {
                $tis_prop_genes{$t}{$g}++;
                $tis_ign_genes{$t}{$g}++; } }

        foreach my $at (keys %{$tis_anc{$t}}) {
            unless(exists $tis_genes{$at}) { next; }
            foreach my $g (keys %{$tis_genes{$at}}) {
                $tis_ign_genes{$t}{$g}++; } }

        unless(-e $tis_exp_dat) {
            my @tgenes = keys %tis_prop_genes;

            for(my $i=0; $i<$#tgenes; $i++) {
                my $g1 = $tgenes[$i];
                if($ispfc) { if(exists $ubiq_genes{$g1}) { next; } }

                for(my $j=($i+1); $j<=$#tgenes; $j++) {
                    my $g2 = $tgenes[$j];
                    if($ispfc) { if(exists $ubiq_genes{$g2}) { next; } }

                    if((exists $ubiq_genes{$g1}) and (exists $ubiq_genes{$g2})) {
                        next; }
                }
            }
        }
    }

    unless(-e $tis_exp_dat) { close TEXP; }
}


# Recording positive edges
open DH, "$ipos" or die "Can't open $ipos!";
chomp(my @dat=<DH>); close DH;

my %within_pose = (); my %tis_pose = ();
my %temp_tis_genes;

$time = runtime(); print "\n$time: Parsing positives ...\n";
foreach my $gpair (@dat) {
    @p = split '\t', $gpair;
    
    # Reasons to skip this edge
    if($p[0] eq $p[1]) { next; }
    unless(exists $all_genes{$p[0]}) { next; }
    unless(exists $all_genes{$p[1]}) { next; }
    if((exists $ubiq_genes{$p[0]}) and (exists $ubiq_genes{$p[1]})) {
        next; }
    if($ispfc) {
        if((exists $ubiq_genes{$p[0]}) or (exists $ubiq_genes{$p[1]})) {
            next; } }

    # Annotating edge to specific tissues
    $e = join '__', sort($p[0], $p[1]);

    if($ionto eq 'genes') {
        %temp_tis_genes = %tis_prop_genes; }
    elsif($ionto eq 'edges') {
        %temp_tis_genes = %tis_genes; }

    foreach my $t (keys %temp_tis_genes) {
        if($ispfc or $inouue) {
            if(((exists $temp_tis_genes{$t}{$p[0]}) and
                    (exists $temp_tis_genes{$t}{$p[1]})) or
                ((exists $temp_tis_genes{$t}{$p[1]}) and
                    (exists $temp_tis_genes{$t}{$p[0]}))) {
                $within_pose{$e}++;
                if($ionto eq 'edges') { $tis_pose{$t}{$e}++; } } }
        else {
            if(((exists $temp_tis_genes{$t}{$p[0]}) and
                    ((exists $temp_tis_genes{$t}{$p[1]}) or
                        (exists $ubiq_genes{$p[1]}))) or
                ((exists $temp_tis_genes{$t}{$p[1]}) and
                    ((exists $temp_tis_genes{$t}{$p[0]}) or
                        (exists $ubiq_genes{$p[0]})))) {
                $within_pose{$e}++;
                if($ionto eq 'edges') { $tis_pose{$t}{$e}++; } } } } }

(my $iposwa = $ipos) =~ s/\.dat/_within\.dat/g; $iposwa =~ s/^.*\///g;
unless(-e $iposwa) {
    open POSW, ">$iposwa";
    foreach my $e (keys %within_pose) {
        @p = split '__', $e;
        print POSW "$p[0]\t$p[1]\t1\n"; }
    close POSW; }

(my $iposg = $ipos) =~ s/\.dat/_within\.genes/g; $iposg =~ s/^.*\///g;
unless(-e $iposg) {
    `Dat2Dab -i $iposwa -E > $iposg`; }


# Recording negative edges within all genes
$time = runtime(); print "\n$time: Recording subsets of negatives ...\n";

(my $inegwa = $ineg) =~ s/\.dab/_within-all\.dab/g; $inegwa =~ s/^.*\///g;
unless (-e $inegwa) {
    `Dat2Dab -i $ineg -g $iposg -X $iubq_glist -o $inegwa`; }

# ... between ubiquitous genes
# ($inegwu = $ineg) =~ s/\.dab/_within-ubiq\.dab/g; $inegwu =~ s/^.*\///g;
# unless (-e $inegwu) {
#     `Dat2Dab -i $inegwa -g $iubq_glist -o $inegwu`; }


# Recording negative edges directly within tissues
(my $inegdt = $ineg) =~ s/\.dab/_dir\.dab/g; $inegwa =~ s/^.*\///g;

# Recording edges for edge-propagation
my %tis_prop_pose = (); my %tis_ign_pose = ();

if($ionto eq 'edges') {
    $time = runtime(); print "\n$time: Propagating tissue edges ...\n";

    # Recording propagated tissue-specific positive edges
    foreach my $t (keys %tis_pose) {
        foreach my $e (keys %{$tis_pose{$t}}) {
            $tis_prop_pose{$t}{$e}++;
            $tis_ign_pose{$t}{$e}++; }

        foreach my $dt (keys %{$tis_des{$t}}) {
            unless(exists $tis_pose{$dt}) { next; }
            foreach my $e (keys %{$tis_pose{$dt}}) {
                $tis_prop_pose{$t}{$e}++;
                $tis_ign_pose{$t}{$e}++; } }

        foreach my $at (keys %{$tis_anc{$t}}) {
            unless(exists $tis_pose{$at}) { next; }
            foreach my $e (keys %{$tis_pose{$at}}) {
                $tis_ign_pose{$t}{$e}++; } } }
    
    # Recording direct tissue-specific negative edges
    $time = runtime(); print "\n$time: Parsing direct negatives ...\n";
    `mkdir -p ./dir-neg`;

    foreach my $t (keys %tis_genes) {
        if(scalar keys %{$tis_genes{$t}} == 0) { next; }

        my $expg = $t.'.exp.genes';
        unless(-e $expg) {
            open EXPG, ">$expg";
            foreach my $g (keys %{$tis_genes{$t}}) {
                print EXPG "$g\n"; }

            unless($ispfc or $inouue) {
                foreach my $g (keys %ubiq_genes) {
                    print EXPG "$g\n"; } }
            close EXPG; }

        my $tc3e = './dir-neg/'.$t.'.dir-neg.dat';
        unless(-e $tc3e) {
            $time = runtime(); print "\t$time: $t\n";
            `Dat2Dab -i $inegwa -g $expg -V 3 -o $tc3e`; } }

    `rm -f tmp.dat`;
    `cat ./dir-neg/*.dir-neg.dat >> tmp.dat`;
    `Dat2Dab -i tmp.dat -o $inegdt`;
}


# Filtering tissues
my %filt_tis = ();
if($ifilt) {
    open FILT, "$ifilt" or die "Can't open $ifilt!";
    while (<FILT>) {
        if($_ =~ /^#/) { next; }
        chomp($_); @p = split '\t', $_;
        $filt_tis{$p[0]} = $p[1]; }
    close FILT; }


# Defining classes: Interaction (I), Tissue-exp (T)
# Class 1: I=1 T=1
# Class 2: I=1 T=0
# Class 3: I=0 T=1
# Class 4: I=0 T=0
$time = runtime(); print "\n$time: Getting tissue-specific gold-standards ...\n";

my ($c1e, $c2e, $c3e, $c4e, $nc1, $nc2, $nc3, $nc4);
my ($ntg, $ntgt, $ntgu, $nc1tt, $nc1tu, $nc2tt, $nc2tu);
my ($nc3tt, $nc3tu, $nc4tt, $nc4tu, $nneg);
my $num_tissues = 0;

print STAT "#Tis.ID\tTis.Desc\tnTisG\tnTisG.t\tnTisG.u\t";
print STAT "nC1\tnC1.tt\tnC1.tu\tnC2\tnC2.tt\tnC2.tu\t";
print STAT "nC3\tnC3.tt\tnC3.tu\tnC3\tnC4.tt\tnC4.tu\tnNeg\n";

open UTS, "$iunp_tis_gmt" or die "Can't open $iunp_tis_gmt!";
chomp(my @tis_gmt = <UTS>); close UTS;

TISSUE: foreach my $gset (@tis_gmt) {
    @p = split '\t', $gset;
    $t = shift(@p); $desc = shift(@p); $desc =~ s/ \([0-9]*\)$//g;
    
    if($t eq $root) { next TISSUE; }
    if($desc =~ /^(hair|milk|saliva|tear|urine|cerebrospinal fluid)$/) { next TISSUE; }
    if($ifilt) { unless(exists $filt_tis{$t}) { next; } }

    $num_tissues++;
    $time = runtime(); print "\t$time:\t$num_tissues\t$t\t$desc\tG:";

    if($ionto eq 'genes') {
        print scalar keys %{$tis_prop_genes{$t}}; }
    if($ionto eq 'edges') {
        print scalar keys %{$tis_genes{$t}}, "\tE:";
        print scalar keys %{$tis_pose{$t}}, "\tpE:";
        print scalar keys %{$tis_prop_pose{$t}}; }
    print "\n";

    $ntg = $ntgt = $ntgu = 0;
    my $expg = $t.'.exp.genes';
    my $uneg = $t.'.une.genes';
    my $igng = $t.'.ign.genes';

    if($ionto eq 'genes') {
        open EXPG, ">$expg";
        foreach my $g (keys %{$tis_prop_genes{$t}}) {
            $ntg++;
            if($inouue) {
                if(exists $ubiq_genes{$g}) { $ntgu++; }
                else { $ntgt++; } }
            else { $ntgt++; }
            print EXPG "$g\n"; }

        open IGNG, ">$igng";
        foreach my $g (keys %{$tis_ign_genes{$t}}) {
            print IGNG "$g\n"; }

        unless($ispfc or $inouue) {
            foreach my $g (keys %ubiq_genes) {
                $ntg++; $ntgu++;
                print EXPG "$g\n";
                print IGNG "$g\n"; } }
        close EXPG;
        close IGNG;

        open UNEG, ">$uneg";
        foreach my $g (keys %all_tspfc_genes) {
            if(exists $tis_ign_genes{$t}{$g}) { next; }
            print UNEG "$g\n"; }
        close UNEG; }

    if($ionto eq 'edges') {
        foreach my $g (keys %{$tis_genes{$t}}) {
            $ntg++;
            if($inouue) {
                if(exists $ubiq_genes{$g}) { $ntgu++; }
                else { $ntgt++; } }
            else { $ntgt++; } }

        unless($ispfc or $inouue) {
            foreach my $g (keys %ubiq_genes) {
                $ntg++; $ntgu++; } } }

    if($ntg == 0) {
        print "\t\tnot enough tissue genes\n";
        `rm -f $t*genes`; next; }

    $nc1 = $nc2 = $nc3 = $nc4 = 0;
    $nc1tt = $nc1tu = $nc2tt = $nc2tu = $nc3tt = $nc3tu = $nc4tt = $nc4tu = 0;

    $c1e = $t.'.c1.dat';
    $c2e = $t.'.c2.dat';
    $c3e = $t.'.c3.dat';
    $c4e = $t.'.c4.dat';

    if($ionto eq 'edges') { 
        # C1
        $time = runtime(); print "\t\t$time: C1 ...\n";
        
        open C1H, ">$c1e";
        foreach my $e (keys %{$tis_prop_pose{$t}}) {
            @r = split '__', $e;
            print C1H "$r[0]\t$r[1]\t1\n";
            $nc1++;

            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc1tu++; }
            else { $nc1tt++; } }
        close C1H;

        # C2
        $time = runtime(); print "\t\t$time: C2 ...\n";

        open C2H, ">$c2e";
        foreach my $e (keys %within_pose) {
            if(exists $tis_ign_pose{$t}{$e}) { next; }

            @r = split '__', $e;
            print C2H "$r[0]\t$r[1]\t2\n";
            $nc2++;

            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc2tu++; }
            else { $nc2tt++; } }
        close C2H;

        # C3
        $time = runtime(); print "\t\t$time: C3 ...\n";

        `cat ./dir-neg/$t.dir-neg.dat > $c3e`;
        foreach my $dt (keys %{$tis_des{$t}}) {
            unless(exists $tis_genes{$dt}) { next; }
            `cat ./dir-neg/$dt.dir-neg.dat >> $c3e`; }
        `Dat2Dab -i $c3e > tmp.dat; mv tmp.dat $c3e`;

        open C3H, "$c3e" or die "Can't open $c3e!";
        while (<C3H>) {
            chomp($_); @r = split '\t', $_;
            $nc3++;
            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc3tu++; }
            else { $nc3tt++; } }
        close C3H;

        # C4
        $time = runtime(); print "\t\t$time: C4 ...\n";

        my $igne = $t.'.ign.dat'; `cp $c3e $igne`;

        foreach my $at (keys %{$tis_anc{$t}}) {
            unless(exists $tis_pose{$at}) { next; }
            `cat ./dir-neg/$at.dir-neg.dat >> $igne`; }
        `Dat2Dab -i $igne > tmp.dat; mv tmp.dat $igne`;
        `Dat2Dab -i $inegdt -x $igne -V 4 -o $c4e`;

        open C4H, "$c4e" or die "Can't open $c4e!";
        while (<C4H>) {
            chomp($_); @r = split '\t', $_;
            $nc4++;
            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc4tu++; }
            else { $nc4tt++; } }
        close C4H;
    }

    elsif($ionto eq 'genes') {
        # C1
        $time = runtime(); print "\t\t$time: C1 ...\n";
        `Dat2Dab -i $iposwa -g $expg -V 1 -o $c1e`;
        
        open C1H, "$c1e" or die "Can't open $c1e!";
        while (<C1H>) {
            chomp($_); @r = split '\t', $_;
            $nc1++;
            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc1tu++; }
            else { $nc1tt++; } }
        close C1H;

        # C2
        $time = runtime(); print "\t\t$time: C2 ...\n";
        `Dat2Dab -i $iposwa -X $igng -V 0 -o $c2e`;
        # `Dat2Dab -i $iposwa -X $igng -V 2 -o $c2e`;

        open C2H, "$c2e" or die "Can't open $c2e!";
        while (<C2H>) {
            chomp($_); @r = split '\t', $_;
            $nc2++;
            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc2tu++; }
            else { $nc2tt++; } }
        close C2H;

        # C3
        $time = runtime(); print "\t\t$time: C3 ...\n";
        `Dat2Dab -i $inegwa -g $expg -V 0 -o $c3e`;
        # `Dat2Dab -i $inegwa -g $expg -V 3 -o $c3e`;

        open C3H, "$c3e" or die "Can't open $c3e!";
        while (<C3H>) {
            chomp($_); @r = split '\t', $_;
            $nc3++;
            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc3tu++; }
            else { $nc3tt++; } }
        close C3H;

        # C4
        $time = runtime(); print "\t\t$time: C4 ...\n";
        `Dat2Dab -i $inegwa -X $igng -V 0 -o $c4e`;
        # `Dat2Dab -i $inegwa -X $igng -V 4 -o $c4e`;

        open C4H, "$c4e" or die "Can't open $c4e!";
        while (<C4H>) {
            chomp($_); @r = split '\t', $_;
            $nc4++;
            if((exists $ubiq_genes{$r[0]}) or (exists $ubiq_genes{$r[1]})) {
                $nc4tu++; }
            else { $nc4tt++; } }
        close C4H;
    }

    @p = sort {$a <=> $b} ($nc2, $nc3, $nc4);
    $nneg = $p[0];

    # Printing gold-std edges
    $out_gsdat = $desc.'.dat';
    `cat $c1e > $out_gsdat`;
    `perl -MList::Util -e 'print List::Util::shuffle <>' $c2e | head -$nneg >> $out_gsdat`;
    `perl -MList::Util -e 'print List::Util::shuffle <>' $c3e | head -$nneg >> $out_gsdat`;
    `perl -MList::Util -e 'print List::Util::shuffle <>' $c4e | head -$nneg >> $out_gsdat`;

    $out_gsdab = $desc.'.dab';
    `Dat2Dab -i $out_gsdat -o $out_gsdab`;

    print STAT "$t\t$desc\t$ntg\t$ntgt\t$ntgu";
    print STAT "\t$nc1\t$nc1tt\t$nc1tu\t$nc2\t$nc2tt\t$nc2tu";
    print STAT "\t$nc3\t$nc3tt\t$nc3tu\t$nc4\t$nc4tt\t$nc4tu\t$nneg\n";

    # last;
    `rm -f $t.c[1-4].dat $t*genes $out_gsdat`;
}

close STAT;

$time = runtime(); print "\n$time: DONE\n\n";


# Subroutines
# Parse OBO to return graph object and a hash of all terms
sub parse_obo {
    my $obofile = shift;
    open OBO, "$obofile" or die "Can't open $obofile!";
        chomp(my @obo=<OBO>); close OBO;
        
    my $dag = Graph->new;

    my $in = 0;
    my ($term, $name, $syn, $rel, $parent, @q);
    LINE: foreach my $line (@obo) {
        if($line =~ /^#/) { next LINE; }
        @q = split(' ', $line);

        if ($line =~ /^$/) { $in = 0; next LINE; }
        elsif ($q[0] eq '[Term]') { $in = 1; }
        elsif ($q[0] eq '[Typedef]') { $in = 0; }
        elsif ($in and ($q[0] eq 'id:')) { ($term = $q[1]) =~ s/:/_/g; }
        elsif ($in and ($q[0] eq 'is_a:')) {
            ($parent = splice(@q, 0, 2)) =~ s/:/_/g;
            $dag->add_edge($term, $parent); }
        elsif ($in and ($q[0] eq 'relationship:')) {
            $rel = $q[1]; ($parent = splice(@q, 0, 3)) =~ s/:/_/g;

            if($rel eq 'has_part') { next LINE; }
            elsif(($rel eq 'part_of') or ($rel =~ /regulates/)) {
                $dag->add_edge($term, $parent); }
            elsif(($rel eq 'develops_from') or ($rel eq 'related_to')) {
                $dag->add_edge($term, $parent); } } }

    return $dag; }



__END__

=head1

Constructs gold-standard across tissues given a template of global gold-std
edges and tissue-specific gene-expression.

=head1 USAGE

construct_tspfc-undir_gold-std.pl [--ipos GOLD-STD_POSITIVES_DAT] [--ineg
GOLD-STD_NEGATIVES_DAT/DAB] [--iunp_tis_gmt TISSUE-GENE_GMT] [--iubq_glist
UBIQ-GENES_LIST] [--ibg_glist BACKGROUND-GENES_LIST] [--ispfc] [--help]

=head1 DESCRIPTION

This script takes in a global template of gold-std edges and tissue-specific
gene-expression and constructs gold-std edges per tissue in 4 classes based on
interaction/tissue combinations.

=head1 ARGUMENTS

=over 12

=item C<--ipos>

Global gold-std positive edges in DAT format.

=item C<--ineg>

Global gold-std negative edges in DAT/DAB format.

=item C<--iunp_tis_gmt>

(Optional) Direct (Unpropagated) tissue-specific gene-expression annotations in
GMT format. Default:
/home/arjunk/data/gene-expression/hprd/hprd-brenda/human_tissue-specific_gene-expression_unprop.gmt

=item C<--iubq_glist>

(Optional) List of ubiquitous genes.
Default: /home/arjunk/data/gene-expression/ubiquitous_genes/human/alldb_ubiquitous-genes_human.txt

=item C<--iobo>

(Optional) Tissue ontology. Default: /home/arjunk/data/ontologies/brenda/BrendaTissue.obo

=item C<--ibg_glist>

(Optional) List of background genes.
Default: /home/arjunk/projects/human-tissue-fln/standards/human_standard.genes

=item C<--ispfc>

(Optional) If provided, consider only tissue-specific genes for edges and ignore
ubiquitous genes for edges.

=item C<--inouue>

(Optional) If provided, consider tissue-gene expression standards as-is, and use
the ubiquitous gene list only to remove edges between ubiquitous genes.

=item C<--ionto>

(Optional) 'genes' or 'edges'.

'genes': Use the ontology to propagate genes based on direct
tissue-gene annotations. This way, genes directly annotated as 'expressed' in a
given tissue or any of its descendants are assigned to that tissue. Also, genes
directly annotated to a given tissue, any of its descendants or any of its
ancestors are avoided from the list of 'unexpressed' genes w.r.t that tissue.

'edges': Use the ontology to propagate edges based on direct
tissue-gene annotations and +ve/–ve standard edges. This way, edges incident on
gene pairs directly annotated as 'expressed' in a given tissue or any of its
descendants are assigned to that tissue. Also, edges incident on gene pairs
directly annotated to a given tissue, any of its descendants or any of its
ancestors are avoided from the list of 'unexpressed' edges w.r.t that tissue.

Default 'genes'.

=item C<--help>

Prints this help message.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut

