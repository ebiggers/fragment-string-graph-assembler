#include "BaseVecVec.h"
#include "Overlap.h"
#include <getopt.h>
#include "Kmer.h"

#include <unordered_map>

class KmerOccurrence {
public:
	unsigned long _read_id	: 32;
	unsigned long _read_pos	: 31;
	unsigned long _rc	: 1;

	KmerOccurrence(unsigned long read_id, unsigned long read_pos, bool rc)
		: _read_id(read_id), _read_pos(read_pos), _rc(rc)
	{ }
};

template <unsigned K>
static void compute_overlaps(const BaseVecVec &bvv, OverlapVecVec &ovv,
			     int overlap_len, int max_edits)
{
	typedef std::unordered_map<Kmer<K>, std::list<KmerOccurrence> > KmerOccurrenceMap;

	ovv.clear();
	ovv.resize(bvv.size());
	KmerOccurrenceMap occ_map;

	for (size_t i = 0; i < bvv.size(); i++) {
		const BaseVec &bv = bvv[i];
		Kmer<K> fwd_kmer;
		Kmer<K> rev_kmer;
		Kmer<K> *kmer;
		bool is_rc;
		if (bv.size() < K)
			continue;
		for (size_t j = 0; j < K - 1; j++) {
			fwd_kmer.push_back(bv[j]);
			rev_kmer.push_rc_front(bv[j]);
		}
		for (size_t j = K; j < bv.size(); j++) {
			fwd_kmer.push_back(bv[j]);
			rev_kmer.push_rc_front(bv[j]);
			if (fwd_kmer < rev_kmer) {
				kmer = &fwd_kmer;
				is_rc = false;
			} else {
				kmer = &rev_kmer;
				is_rc = true;
			}
			occ_map[*kmer].push_back(KmerOccurrence(i, j, is_rc));
		}
	}
}

static const char *optstring = "l:e:h";
static const struct option longopts[] = {
	{"overlap-len", required_argument, NULL, 'l'},
	{"max-edits",   required_argument, NULL, 'e'},
	END_LONGOPTS
};

DEFINE_USAGE(
"Usage: compute-overlaps BVV_FILE OVERLAPS_FILE\n"
"\n"
"  -l, --overlap-len=LEN\n"
"  -e, --max-edits=MAX_EDITS\n"
"  -h, --help\n"
);

int main(int argc, char *argv[])
{
	int c;
	int overlap_len = 25;
	int max_edits = 0;
	for_opt(c) {
		switch (c) {
		case 'l':
			overlap_len = parse_long(optarg, "--overlap-len");
			break;
		case 'e':
			max_edits = parse_long(optarg, "--max-edits");
			break;
		PROCESS_OTHER_OPTS
		}
	}
	argc -= optind;
	argv += optind;
	USAGE_IF(argc != 2);

	if (max_edits != 0)
		unimplemented();

	BaseVecVec bvv(argv[1]);
	OverlapVecVec ovv;
	compute_overlaps<24>(bvv, ovv, overlap_len, max_edits);
	ovv.write(argv[2]);
}
