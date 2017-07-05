open FH, "$ARGV[0]"; # gene_association_yeast.mf.closed
open GH, "/Genomics/Users/arjunk/data/functional-annotations/go/gene_ontology.all.terms";
($fo = $ARGV[0]) =~ s/gene_association_([a-z\-]*)(\.*[A-Z\-]*)\.(..)\.(dir|closed)/go\3_\1\2.\4.gmt/g;
print $fo, "\n";
open HH, ">$fo";

%gs_desc = ();
while(<GH>) {
    chomp($_); @p = split '\t', $_;
    $gs_desc{$p[0]} = $p[1]; }

%gs_genes = ();
while(<FH>) {
    chomp($_); @p = split '\t', $_;
    $gs_genes{$p[0]}{$p[1]}++; }

%gs_size = ();
foreach (keys %gs_genes) {
    $gs_size{$_} = scalar keys %{$gs_genes{$_}}; }

foreach $gs (sort {$gs_size{$b} <=> $gs_size{$a}} keys %gs_genes) {
    print HH "$gs\t$gs_desc{$gs} (", scalar keys %{$gs_genes{$gs}},")";
    
    foreach $gene (keys %{$gs_genes{$gs}}) {
        print HH "\t$gene"; }
    print HH "\n"; }

close HH;
