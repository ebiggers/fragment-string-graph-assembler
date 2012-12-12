#!/usr/bin/env perl

use strict;
use Getopt::Long;
$| = 1;

my $help = 0;
my $nbytes = -1;
my $nlines = -1;

my $USAGE = "Usage: fasta_head.pl [--bytes BYTES] [--lines LINES]\n";

my $res = GetOptions("help"        => \$help,
                     "bytes=n"     => \$nbytes,
                     "lines=n"     => \$nlines);

if (!$res) {
    die "$USAGE";
}

if ($help) {
    print "$USAGE";
    exit 0;
}

if ($nbytes == -1 && $nlines == -1) {
    $nlines = 10;
}

my $bytes = 0;
my $lines = 0;

while (<>) {
    chomp;
    if (/^>/) {
        print "$_\n";
    } else {
        if ($lines++ == $nlines) {
            last;
        }
        if ($nbytes != -1 && $bytes + length($_) >= $nbytes) {
            my $s = substr($_, 0, $nbytes - $bytes);
            print "$s\n";
            last;
        }
        print "$_\n";
        $bytes += length($_);
    }
}
