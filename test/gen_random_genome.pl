#!/usr/bin/env perl

my $USAGE =
"Usage: gen_random_genome.pl LENGTH\n";

my $length = $ARGV[0] or die $USAGE;

print ">random_genome\n";

my $j = 0;
srand(1);
my %bases = (
    0 => "A",
    1 => "C",
    2 => "G",
    3 => "T"
);
for (my $i = 0; $i < $length; $i++) {
    my $n = int(rand(4));
    my $base = $bases{$n};
    print $base;
    if (++$j == 70) {
        $j = 0;
        print "\n";
    }
}
