#!/usr/bin/perl

use strict;
use Getopt::Long;
$| = 1;

my $help = 0;
my $read_len = 100;
my $read_sep = 10;


my $USAGE = "Usage: simulate_uniform_reads.pl [--help] [--read-len=LEN]
[--read-sep=SEP] GENOME_FASTA\n";

my $res = GetOptions("help"        => \$help,
                     "read-len=n", => \$read_len,
                     "read-sep=n", => \$read_sep);

if (!$res) {
    die "$USAGE";
}

if ($help) {
    print "$USAGE";
    exit 0;
}
                     


my $file = $ARGV[0];
open FILE, "<", $file or die "Cannot open \"$file\": $!";

my $seq = "";
while (<FILE>) {
    chomp;
    if (/^>/) {
        next;
    } else {
        $seq .= $_;
    }
}
print "$seq\n";

for (my $i = 0; ; $i += $read_sep) {
    if ($i + $read_len > length($seq)) {
        last;
    }
    my $read = substr($seq, $i, $read_len);
    print ">\n$read\n";
}

