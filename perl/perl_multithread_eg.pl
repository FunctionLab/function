#!/opt/local/bin/perl -w
use threads;
use strict;
use warnings;

my @a = ();
my @b = ();

sub sleeping_sub ( $ $ $ );

print "Starting main program\n";

my $nb_process = 10;
my $nb_compute = 20;
my $i=0;
my @running = ();
my @Threads;
while (scalar @Threads < $nb_compute) {
 	@running = threads->list(threads::running);
	print "LOOP $i\n";
	print "  - BEGIN LOOP >> NB running threads = ".(scalar @running)."\n";

	if (scalar @running < $nb_process) {
 		my $thread = threads->new( sub { sleeping_sub($i, \@a, \@b) });
		push (@Threads, $thread);
		my $tid = $thread->tid;
		print "  - starting thread $tid\n";
	}
	@running = threads->list(threads::running);
	print "  - AFTER STARTING >> NB running Threads = ".(scalar @running)."\n";
	foreach my $thr (@Threads) {
		if ($thr->is_running()) {
			my $tid = $thr->tid;
			print "  - Thread $tid running\n";
		}
		elsif ($thr->is_joinable()) {
			my $tid = $thr->tid;
			$thr->join;
			print "  - Results for thread $tid:\n";
			print "  - Thread $tid has been joined\n";
		}
	}

	@running = threads->list(threads::running);
	print "  - END LOOP >> NB Threads = ".(scalar @running)."\n";
	$i++;
}
print "\nJOINING pending threads\n";
while (scalar @running != 0) {
	foreach my $thr (@Threads) {
		$thr->join if ($thr->is_joinable());
	}
	@running = threads->list(threads::running);
}

print "NB started threads = ".(scalar @Threads)."\n";
print "End of main program\n";

sub sleeping_sub ( $ $ $ ) {
	sleep(4);
}
