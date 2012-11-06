#!/usr/bin/perl

use strict;
use Getopt::Long;
$| = 1;

my $help = 0;
my $read_len = 100;
my $read_sep = 10;
my $allow_rc = 0;


my $USAGE =
"Usage: simulate_uniform_reads.pl [--help] [--read-len=LEN]
        [--read-sep=SEP] [--allow-rc] GENOME_FASTA\n";

my $res = GetOptions("help"        => \$help,
                     "read-len=n", => \$read_len,
                     "read-sep=n", => \$read_sep,
                     "allow-rc",   => \$allow_rc);

if (!$res) {
    die "$USAGE";
}

if ($help) {
    print "$USAGE";
    exit 0;
}

my $file = $ARGV[0] or die $USAGE;

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

my $read_no = 1;
for (my $i = 0; ; $i += $read_sep) {
    if ($i + $read_len > length($seq)) {
        last;
    }
    my $read = substr($seq, $i, $read_len);
    $read =~ tr/acgt/ACGT/;
    if ($allow_rc && int(rand(2)) == 0) {
        $read = reverse $read;
        $read =~ tr/ACGT/TGCA/;
    }
    print ">read_$read_no\n$read\n";
    $read_no++;
}

