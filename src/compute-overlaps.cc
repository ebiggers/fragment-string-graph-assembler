#include "BaseVecVec.h"
#include "Overlap.h"
#include <getopt.h>
#include "Kmer.h"

#include <unordered_map>

class KmerOccurrence {
private:
	unsigned long _read_id	: 32;
	unsigned long _read_pos	: 31;
	unsigned long _rc	: 1;
public:

	KmerOccurrence(unsigned long read_id, unsigned long read_pos, bool rc)
		: _read_id(read_id), _read_pos(read_pos), _rc(rc)
	{ }

	unsigned long get_read_id() const {
		return _read_id;
	}

	unsigned long get_read_pos() const {
		return _read_pos;
	}

	bool is_rc() const {
		return _rc;
	}
};

static void assert_seed_valid(const BaseVec & bv1, const BaseVec & bv2,
			      const unsigned pos1, const unsigned pos2,
			      const unsigned len,
			      const bool is_rc1, const bool is_rc2)
{
	unsigned i;
	if (is_rc1 == is_rc2) {
		for (i = 0; i < len; i++)
			if (bv1[pos1 + i] != bv2[pos2 + i])
				goto seed_invalid;
	} else if (is_rc2) {
		for (i = 0; i < len; i++)
			if (bv1[pos1 + i] != (3 ^ bv2[pos2 + len - 1 - i]))
				goto seed_invalid;
	} else if (is_rc1) {
		for (i = 0; i < len; i++)
			if (bv1[pos1 + len - 1 - i] != (3 ^ bv2[pos2 + i]))
				goto seed_invalid;
	} else {
		unreachable();
	}
	return;
seed_invalid:
	fatal_error("SEED INVALID (pos1 = %u, pos2 = %u, len = %u, "
		    "is_rc1 = %d, is_rc2 = %d)",
		    pos1, pos2, len, is_rc1, is_rc2);
}

static void extend_match(const BaseVec & bv1,
			 const BaseVec & bv2,
			 unsigned & pos1,
			 unsigned & pos2,
			 unsigned & len,
			 const bool is_rc1,
			 const bool is_rc2)
{
	if (!is_rc1 && !is_rc2) {
	}
}

static bool find_overlap(const BaseVecVec & bvv,
			 const KmerOccurrence occ1,
			 const KmerOccurrence occ2,
			 const unsigned min_overlap_len,
			 const unsigned max_edits,
			 const unsigned K,
			 Overlap &o)
{
	/* ------------>
	 *         ------------->
	 *
	 *
	 *  -------------->
	 *           <---------------
	 *
	 *
	 *           -------------->
	 *  ------------->
	 *
	 *
	 *          --------------->
	 *  <-------------
	 */
	const BaseVec & bv1 = bvv[occ1.get_read_id()];
	const BaseVec & bv2 = bvv[occ2.get_read_id()];
	unsigned pos1 = occ1.get_read_pos();
	unsigned pos2 = occ2.get_read_pos();
	unsigned len = K;
	const bool is_rc1 = occ1.is_rc();
	const bool is_rc2 = occ2.is_rc();
	assert_seed_valid(bv1, bv1, pos1, pos2, len, is_rc1, is_rc2);
	extend_match(bv1, bv2, pos1, pos2, len, is_rc1, is_rc2);
	if (len >= min_overlap_len) {
		unsigned long read_1_beg = pos1;
		unsigned long read_1_end = pos1 + len - 1;
		unsigned long read_2_beg = pos2;
		unsigned long read_2_end = pos2 + len - 1;

		unsigned num_extremes = 0;
		if (read_1_beg == 0)
			num_extremes++;
		if (read_1_end == bv1.size() - 1)
			num_extremes++;
		if (read_2_beg == 0)
			num_extremes++;
		if (read_2_end == bv2.size() - 1)
			num_extremes++;

		if (num_extremes >= 2) {
			if (is_rc1)
				std::swap(read_1_beg, read_1_end);
			if (is_rc2)
				std::swap(read_2_beg, read_2_end);

			o.set(occ1.get_read_id(), read_1_beg, read_1_end,
			      occ2.get_read_id(), read_2_beg, read_2_end);
			return true;
		}
	}
	return false;
}

template <unsigned K>
static void load_kmer_occurrences(const BaseVecVec &bvv,
				  std::unordered_map<Kmer<K>,
				  		     std::vector<KmerOccurrence> > &occ_map)
{
	occ_map.clear();
	unsigned long num_kmer_occurrences = 0;
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
			rev_kmer.push_front(bv[j] ^ 3);
		}
		for (size_t j = K - 1; j < bv.size(); j++) {
			fwd_kmer.push_back(bv[j]);
			rev_kmer.push_front(bv[j] ^ 3);
			if (fwd_kmer < rev_kmer) {
				kmer = &fwd_kmer;
				is_rc = false;
			} else {
				kmer = &rev_kmer;
				is_rc = true;
			}
			occ_map[*kmer].push_back(KmerOccurrence(i, j - (K - 1), is_rc));
			num_kmer_occurrences++;
		}
	}
	info("Loaded %lu %lu-mer occurrences into hash map",
	     num_kmer_occurrences, K);
}


template <unsigned K>
static void compute_overlaps(const BaseVecVec &bvv, 
			     const unsigned min_overlap_len,
			     const unsigned max_edits,
			     OverlapVecVec &ovv)
{
	typedef std::vector<KmerOccurrence> KmerOccurrenceVec;
	typedef std::unordered_map<Kmer<K>, KmerOccurrenceVec > KmerOccurrenceMap;

	ovv.clear();
	ovv.resize(bvv.size());

	KmerOccurrenceMap occ_map;

	load_kmer_occurrences(bvv, occ_map);

	Overlap o;
	unsigned long num_overlaps = 0;
	for (typename KmerOccurrenceMap::const_iterator it = occ_map.begin();
	     it != occ_map.end(); it++)
	{
		const KmerOccurrenceVec & occs = it->second;
		for (size_t i = 0; i < occs.size(); i++) {
			for (size_t j = i + 1; j < occs.size(); j++) {
				if (find_overlap(bvv, occs[i], occs[j],
						 min_overlap_len, max_edits, K, o)
				    && std::find(ovv[i].begin(), ovv[i].end(), o)
				      == ovv[i].end())
				{
					num_overlaps++;
					ovv[i].insert(o);
				}
			}
		}
	}
	info("Found %lu overlaps", num_overlaps);
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
	unsigned overlap_len = 25;
	unsigned max_edits = 0;
	for_opt(c) {
		switch (c) {
		case 'l':
			overlap_len = parse_long(optarg, "--overlap-len",
						 1, UINT_MAX);
			break;
		case 'e':
			max_edits = parse_long(optarg, "--max-edits",
					       0, UINT_MAX);
			break;
		PROCESS_OTHER_OPTS
		}
	}
	argc -= optind;
	argv += optind;
	USAGE_IF(argc != 2);

	if (max_edits != 0)
		unimplemented();

	BaseVecVec bvv(argv[0]);
	OverlapVecVec ovv;
	compute_overlaps<24>(bvv, overlap_len, max_edits, ovv);
	ovv.write(argv[1]);
}
