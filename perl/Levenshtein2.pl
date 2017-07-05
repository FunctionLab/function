# From http://www.mgilleland.com/ld/ldperl2.htm

my ($s1, $s2) = (@ARGV);
print "The Levenshtein distance between $s1 and $s2 is: " . levenshtein($s1, $s2) . "\n";

sub levenshtein($$){
  my @A=split //, lc shift;
  my @B=split //, lc shift;
  my @W=(0..@B);
  my ($i, $j, $cur, $next);
  for $i (0..$#A){
	$cur=$i+1;
	for $j (0..$#B){
		$next=min(
			$W[$j+1]+1,
			$cur+1,
			($A[$i] ne $B[$j])+$W[$j]
		);
		$W[$j]=$cur;
		$cur=$next;
	}
	$W[@B]=$next;
  }
  return $next;
}

sub min($$$){
  if ($_[0] < $_[2]) { pop @_; } else { shift @_; }
  return $_[0] < $_[1]? $_[0]:$_[1];
}

