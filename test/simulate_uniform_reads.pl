#!/usr/bin/env perl

use strict;
use Getopt::Long;
$| = 1;

my $help = 0;
my $read_len = 100;
my $read_sep = 10;
my $allow_rc = 0;
my $random_sep = 0;
my $read_pos_log = undef;
my $coverage = 0;
my $seed = 1;

my $USAGE =
"Usage: simulate_uniform_reads.pl [--help] [--read-len=LEN]
        [--read-sep=SEP] [--allow-rc] [--coverage=X]
        [--read-pos-log=FILE] GENOME_FASTA\n";

my $res = GetOptions("help"           => \$help,
                     "read-len=n"     => \$read_len,
                     "read-sep=n"     => \$read_sep,
                     "allow-rc"       => \$allow_rc,
                     "coverage=f"     => \$coverage,
                     "read-pos-log=s" => \$read_pos_log,
                     "seed=n"         => \$seed);

if (!$res) {
    die "$USAGE";
}

if ($help) {
    print "$USAGE";
    exit 0;
}

srand($seed);

my $file = $ARGV[0] or die $USAGE;

open FILE, "<", $file or die "Cannot open \"$file\": $!";

if ($read_pos_log) {
    open READ_POS_LOG, ">", $read_pos_log or die "Cannot open \"$read_pos_log\": $!";
}

my $seq = "";
while (<FILE>) {
    chomp;
    if (/^>/) {
        next;
    } else {
        $seq .= $_;
    }
}
$seq =~ tr/acgt/ACGT/;

my $read_no = 1;

sub sim_read_at_pos {
    my $pos = shift;

    my $read = substr($seq, $pos, $read_len);
    my $dir;

    if ($allow_rc && int(rand(2)) == 0) {
        $read = reverse $read;
        $read =~ tr/ACGT/TGCA/;
        $dir = "-";
    } else {
        $dir = "+";
    }
    print ">read_$read_no\n$read\n";
    if ($read_pos_log) {
        print READ_POS_LOG "read_$read_no\t$pos\t$dir\n";
    }
    $read_no++;
}

if ($coverage) {
    my $num_reads = length($seq) * $coverage / $read_len;
    for (my $i = 0; $i < $num_reads; $i++) {
        my $pos = int(rand(length($seq) - $read_len));
        if ($pos >= 0) {
            sim_read_at_pos($pos);
        }
    }
} else {
    for (my $pos = 0; $pos + $read_len <= length($seq); $pos += $read_sep) {
        sim_read_at_pos($pos);
    }
}

if ($read_pos_log) {
    close READ_POS_LOG or die "Cannot close \"$read_pos_log\": $!";
}
