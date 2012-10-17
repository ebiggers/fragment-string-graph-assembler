#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <vector>
#include <limits>
#include "Read.h"
#include "dna-util.h"
#include "util.h"
#include "parse-reads.h"
using std::vector;

static const char *optstring = "hO:e:t:";
static const struct option longopts[] = {
	{"overlap-len", required_argument, NULL, 'O'},
	{"errors-per-overlap", required_argument, NULL, 'e'},
	{"threads", required_argument, NULL, 't'},
	{"help", no_argument, NULL, 'h'},
	{NULL, 0, NULL, 0},
};

static void usage()
{
	const char *usage_str = 
"USAGE: assemble (FASTA_FILE | FASTQ_FILE)...";
	puts(usage_str);
}

static void load_read(Read &read, void *vec_reads)
{
	vector<Read> *reads = static_cast<vector<Read>*>(vec_reads);
	reads->push_back(read);
}

struct Overlap {
	unsigned read_1_idx : 24;
	unsigned read_2_idx : 24;
	unsigned read_1_beg : 12;
	unsigned read_1_end : 12;
	unsigned read_2_beg : 12;
	unsigned read_2_end : 12;
};

static void add_overlaps(const vector<Read> &reads, size_t read_1_idx,
			 size_t read_2_idx, size_t overlap_len,
			 size_t errors_per_overlap, vector<Overlap> &overlaps)
{
	const Read &read_1 = reads[read_1_idx];
	const Read &read_2 = reads[read_2_idx];
	const char *seq_1  = read_1.seq;
	size_t seq_1_len   = read_1.seq_len;
	const char *seq_2  = read_2.seq;
	size_t seq_2_len   = read_2.seq_len;
	Overlap overlap;
	overlap.read_1_idx = read_1_idx;
	overlap.read_2_idx = read_2_idx;

	size_t num_differences;
	if (seq_1_len <= overlap_len || seq_2_len <= overlap_len)
		return;

	/*    read_1
	 * -------------->
	 *          ---------------->
	 *                 read_2
	 */
	num_differences = 0;
	for (size_t i = 0; i < overlap_len; i++)
		if (seq_1[seq_1_len - overlap_len + i] != seq_2[i])
			if (++num_differences > errors_per_overlap)
				goto overlap_2;
	overlap.read_1_beg = seq_1_len - overlap_len;
	overlap.read_1_end = seq_1_len - 1;
	overlap.read_2_beg = 0;
	overlap.read_2_end = overlap_len - 1;
	overlaps.push_back(overlap);
overlap_2:
	/*
	 *     read_1
	 * -------------->
	 *          <----------------
	 *                 read_2
	 */
	num_differences = 0;
	for (size_t i = 0; i < overlap_len; i++)
		if (base_ascii_to_bin(seq_1[seq_1_len - overlap_len + i])
		    != complement(base_ascii_to_bin(seq_2[seq_2_len - 1 - i])))
			if (++num_differences > errors_per_overlap)
				goto overlap_3;
	overlap.read_1_beg = seq_1_len - overlap_len;
	overlap.read_1_end = seq_1_len - 1;
	overlap.read_2_beg = seq_2_len - 1;
	overlap.read_2_end = seq_2_len - overlap_len;
	overlaps.push_back(overlap);
overlap_3:
	/*
	 *     read_1
	 * <------------- 
	 *          <----------------
	 *                 read_2
	 */
	num_differences = 0;
	for (size_t i = 0; i < overlap_len; i++)
		if (seq_1[overlap_len - 1 - i] != seq_2[seq_2_len - 1 - i])
			if (++num_differences > errors_per_overlap)
				goto overlap_4;
	overlap.read_1_beg = overlap_len - 1;
	overlap.read_1_end = 0;
	overlap.read_2_beg = seq_2_len - 1;
	overlap.read_2_end = seq_2_len - overlap_len;
	overlaps.push_back(overlap);
overlap_4:
	/*
	 *     read_1
	 * <------------- 
	 *          ---------------->
	 *                 read_2
	 */
	num_differences = 0;
	for (size_t i = 0; i < overlap_len; i++)
		if (seq_1[i] != ascii_complement(seq_2[overlap_len - 1 - i]))
			if (++num_differences > errors_per_overlap)
				return;
	overlap.read_1_beg = overlap_len - 1;
	overlap.read_1_end = 0;
	overlap.read_2_beg = 0;
	overlap.read_2_end = overlap_len - 1;
	overlaps.push_back(overlap);
}

static void compute_overlaps(const vector<Read> &reads, size_t overlap_len,
			     size_t errors_per_overlap, vector<Overlap> &overlaps,
			     int num_threads)
{
	size_t num_reads = reads.size();
	overlaps.clear();
	#pragma omp parallel num_threads(num_threads)
	{
		vector<Overlap> my_overlaps;
		#pragma omp for schedule(dynamic, 32)
		for (size_t i = 0; i < num_reads; i++) {
			for (size_t j = 0; j < num_reads; j++) {
				if (i != j) {
					add_overlaps(reads, i, j, overlap_len,
						     errors_per_overlap, my_overlaps);
				}
			}
		}
		#pragma omp critical
		{
			overlaps.insert(overlaps.end(), my_overlaps.begin(),
					my_overlaps.end());
		}
	}
	info("Found %zu overlaps", overlaps.size());
}

int main(int argc, char **argv)
{
	int c;
	size_t overlap_len = 50;
	size_t errors_per_overlap = 2;
	int num_threads = 0;

	while ((c = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
		switch (c) {
		case 'O':
			overlap_len = strtoul(optarg, NULL, 10);
			break;
		case 'e':
			errors_per_overlap = strtoul(optarg, NULL, 10);
			break;
		case 't':
			num_threads = atoi(optarg);
			break;
		case 'h':
		default:
			usage();
			exit((c == 'h') ? 0 : 2);
		}
	}

	argc -= optind;
	argv += optind;
	if (argc == 0) {
		usage();
		return 2;
	}
	if (num_threads == 0)
		num_threads = get_default_num_threads();

	const char **reads_files = const_cast<const char **>(argv);
	int num_reads_files = argc;
	vector<Read> reads;

	info("Launching assembly");
	info("Parameters:");
	info("  overlap_len = %zu", overlap_len);
	info("  errors_per_overlap = %zu", errors_per_overlap);
	info("  reads files:");
	for (int i = 0; i < num_reads_files; i++)
		info("    \"%s\"", reads_files[i]);

	for_all_reads(reads_files, num_reads_files, load_read, &reads, 1);

#if 0
	for (size_t i = 0; i < reads.size(); i++) {
		const Read &r = reads[i];
		for (size_t j = 0; j < r.seq_len; j++)
			putchar(r.seq[j]);
		putchar('\n');
	}
#endif
	ssize_t min_len = std::numeric_limits<ssize_t>::max();
	ssize_t max_len = std::numeric_limits<ssize_t>::min();
	size_t total_len = 0;
	size_t num_reads = reads.size();
	vector<Read>::const_iterator it;
	for (it = reads.begin(); it != reads.end(); it++) {
		min_len = std::min(it->seq_len, min_len);
		max_len = std::max(it->seq_len, max_len);
		total_len += it->seq_len;
	}
	double avg_len = double(total_len) / double(num_reads);

	info("Loaded %zu reads", reads.size());
	info("  [min_len = %zu, avg_len = %g, max_len = %zu]",
	     min_len, avg_len, max_len);

	vector<Overlap> overlaps;
	compute_overlaps(reads, overlap_len, errors_per_overlap,
			 overlaps, num_threads);
	return 0;
}
