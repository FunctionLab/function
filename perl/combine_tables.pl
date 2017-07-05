#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, @itab, $itab_list, $otab);
pod2usage( -exitstatus => 2, -verbose => 2 ) if (@ARGV==0);
GetOptions('help' => \$help,
     'itablist=s' => \$itab_list,
      'itab=s{,}' => \@itab,
         'otab=s' => \$otab) or pod2usage( -exitstatus => 2, -verbose => 2 );
pod2usage( -exitstatus => 2, -verbose => 2 ) if ($help);

open HH, ">$out_file";

my (@p, @q, $tmp_head, $header);
if(scalar(@itab) > 0) {
    foreach (@itab) {
        open FH, "$_" or die "Can't find $_!";
            chomp(my @table = <FH>); close FH;
        
        $header = shift(@table);
        # foreach (@table) {
        # }
    }
}

elsif($itab_list) {
    open GH, "$itab_list" or die "Can't find $itab_list!";
        chomp(my @list = <GH>); close GH;

    $header = '';
    foreach (@list) {
        @p = split '\t', $_;

        open FH, "$p[0]" or die "Can't find $p[0]!";
            chomp(my @table = <FH>); close FH;

        $tmp_head = shift(@table);
        @q = split '\t', $tmp_head;
        @q = map {$_ = $p[1].'.'.$_ } @q;
    }

}

close HH;

__END__

=head1

<Brief_Desc>

=head1 USAGE

./code.pl [--i INPUT_FILE] [--p OPTION] [--o OUTPUT_FILE] [--helpp]

=head1 DESCRIPTION

This script ...

=head1 ARGUMENTS

=over 12

=item C<--i>

Input file

=item C<--p>

Option

=item C<--o>

Output file

=item C<--helpp>

Prints this documentation.

=back

=head1 AUTHOR

Arjun Krishnan <arjunk@princeton.edu>

=head1 DATE

2012 May 15

=cut


$fi=$ARGV[0];		# List of tables to combine
chomp($fo=$ARGV[1]);

open GH,"$fi";
chomp(@files=<GH>);
close GH;

open HH,">$fo";

%allgenes=();
for($i=0;$i<=$#files;$i++)
{
	open FH,"$files[$i]";
	chomp(@f=<FH>);
	close FH;

	$numgenes=0;
	foreach(@f)
	{
		$_=~s///g;

		@p=();
		@p=split("\t",$_);

		$allgenes{$p[0]}++;
		$numgenes++;
	}
	
	print "$files[$i]\t$numgenes\n";
}

print HH "Gene";
%geneattribute=();
for($i=0;$i<=$#files;$i++)
{
	open FH,"$files[$i]";
	chomp(@f=<FH>);
	close FH;

	$col=''; $col=$files[$i];
	$col=~s/\.txt$//g; $col=~s/\.genesets$//g;
	print HH "\t$col";

	%tempgeneatt=();
	foreach(@f)
	{
		$_=~s///g;

		@p=();
		@p=split("\t",$_);

		$tempgeneatt{$p[0]}=$p[1];
	}

	foreach $gene (keys %allgenes)
	{
		$att='';
		if(exists $tempgeneatt{$gene}) { $att=$tempgeneatt{$gene}; }

		push(@{$geneattribute{$gene}}, $att);
	}
}
print HH "\n";

foreach $gene (sort {$a cmp $b} keys %geneattribute)
{
	print HH "$gene";

	for($j=0;$j<=$#{$geneattribute{$gene}};$j++)
	{
		print HH "\t${$geneattribute{$gene}}[$j]";
	}

	print HH "\n";
}

close HH;
