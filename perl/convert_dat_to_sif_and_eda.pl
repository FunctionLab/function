#################################################################################
# This script converts a network in DAT format - <nodeA> <nodeB> <weight> -	#
# into SIF format and an accompanying edge attribute ATT file, both for		#
# visualizing the network in cytoscape.						#
#################################################################################

# DAT file
# List of nodes (to identify singletons not in the DAT file and add them to the end of the SIF file)
# Edge type tag; will be used in the output files
# Edges with scores above this cutoff will be retained

my ($idat, $ietag, $icut, $igenes) = @ARGV;

(my $osif = $idat) =~ s/\.dat$/\.sif/g;
(my $oattrs = $idat) =~ s/\.dat$/\.attrs/g;
#(my $edatag = $idat) =~ s/\.dat//g;

open SIF,">$osif";
open ATT,">$oattrs";

#print ATT "$edatag\n";
my %net_nodes=();

open DAT, "$idat" or die "Can't open $idat!";
while (<DAT>) {
    if($_ =~ /^#/) { next; }
    chomp($_); my @p = split '\t', $_;

    if($p[2] < $icut) { next; }

    $net_nodes{$p[0]}++; $net_nodes{$p[1]}++;
    print SIF "$p[0]\t$ietag\t$p[1]\n";
    print ATT "$p[0] ($ietag) $p[1] = $p[2]\n"; }
close DAT;

if($igenes) {
    open GEN, "$igenes" or die "Can't open $igenes!";
    while (<GEN>) {
        if($_ =~ /^#/) { next; }
        chomp($_); my @p = split '\t', $_;
        if(exists $net_nodes{$p[0]}) { next; }
        print SIF "$p[0]\n"; }
    close GEN; }

close SIF;
close ATT;
