#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Long;
#use Data::Dumper;
#use Time::SoFar qw(runtime);

if($#ARGV<2) {
	print "\nThis script takes an obo file and outputs all pairs of GO terms that have a 'regulatory' relationsip.";
	print "\nUsage: ./get_go-regulatory-relationships.pl <obo_file> <outfile_reg-target_goterm_pairs> <outfile_list_of_concerned_goterms>\n\n";
	exit;
}

my ($obofile, $outfile1, $outfile2) = @ARGV;
open FH, "$obofile" or die "Can't open file: $obofile!"; chomp(my @obo=<FH>); close FH;
open RH, ">$outfile1";	# This file will contain all the regulatory-target GO term pairs
open RC, ">$outfile2";	# This file will contain a unique set of GO terms from the previous file, with 'regulatory' or 'core' annotation

my($token1, $token2, $id, %term, %namespace, %alt_id, %relation_isa, %relation_partof, %relation_regulates, %function_att); my $numids = 0;
for(my $i=0;$i<=$#obo;$i++) {
	if($obo[$i]=~/^\[Term\]$/) {
		$numids++; # if(($numids%1000)==0) { $time = runtime(); print "$time\t$numids\t$obo[$i+1]\n"; }
		do {
			$i++; $token1 = (split(" ",$obo[$i]))[1];
			if($obo[$i]=~/^id: /) { $id = $token1; }
			elsif($obo[$i]=~/^name: /) { $term{$id} = (split(": ",$obo[$i]))[1]; }
			elsif($obo[$i]=~/^namespace: /) { $namespace{$id} = $token1; }
			elsif($obo[$i]=~/^alt_id: /) { push(@{$alt_id{$id}}, $token1); }
			elsif($obo[$i]=~/^is_obsolete: /) { if($token1 eq 'true') { delete $term{$id}; } }
			
			elsif($obo[$i]=~/^is_a: /) { push(@{$relation_isa{$id}}, $token1); }
			elsif($obo[$i]=~/^relationship: /) {
				$token2 = (split(" ",$obo[$i]))[2];
				if($token1 eq 'part_of') { push(@{$relation_partof{$id}}, $token2); }
				elsif(($token1 eq 'regulates')||($token1 eq 'positively_regulates')||($token1 eq 'negatively_regulates')) {
					push(@{$relation_regulates{$id}}, $token2);
					$function_att{$id} = 'R'; unless(exists $function_att{$token2}) { $function_att{$token2} = 'C'; }
				}
			}
		}until($obo[$i+1]=~/^$/);
	}
}

my %uniq_functions;
foreach my $id (keys %relation_regulates) {
	if(exists $term{$id}) {
		$uniq_functions{$id}++;
		foreach my $target (@{$relation_regulates{$id}}) {
			if(exists $term{$target}) {
				print RH "$function_att{$id}-$function_att{$target}\t$id\t$term{$id}\t$target\t$term{$target}\n";
				$uniq_functions{$target}++;
			}
		}
	}
}

foreach my $id (sort {$function_att{$b} cmp $function_att{$a}} keys %uniq_functions) {
	if(exists $term{$id}) {
		print RC "$function_att{$id}\t$id\t$term{$id}\n";
	}
}

close RH;
close RC;
