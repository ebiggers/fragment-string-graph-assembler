#include "BaseVecVec.h"
#include "Overlap.h"
#include <getopt.h>
#include "Kmer.h"

#include <unordered_map>
#include <ostream>

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
	friend std::ostream & operator<<(std::ostream & os, const KmerOccurrence & occ)
	{
		return os << "KmerOccurrence { _read_id: " << occ._read_id << 
			", _read_pos: " << occ._read_pos <<
			", _rc: " <<occ._rc << "}";
	}
};

static void assert_seed_valid(const BaseVec & bv1,
			      const BaseVec & bv2,
			      const unsigned pos1,
			      const unsigned pos2,
			      const unsigned len,
			      const bool is_rc1, const bool is_rc2,
			      const char *description = "SEED")
{
	unsigned i;
	if (pos1 + len > bv1.size())
		goto seed_invalid;
	if (pos2 + len > bv2.size())
		goto seed_invalid;

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
	}
	return;
seed_invalid:
	std::cerr << bv1 << std::endl;
	std::cerr << bv2 << std::endl;
	fatal_error("%s INVALID (pos1 = %u, pos2 = %u, len = %u, "
		    "is_rc1 = %d, is_rc2 = %d)", description,
		    pos1, pos2, len, is_rc1, is_rc2);
}

static void assert_overlap_valid(const Overlap & o, const BaseVecVec & bvv,
				 const unsigned min_overlap_len,
				 const unsigned max_edits)
{
	if (max_edits > 0)
		unimplemented();
	unsigned long read_1_idx, read_1_beg, read_1_end;
	unsigned long read_2_idx, read_2_beg, read_2_end;
	o.get(read_1_idx, read_1_beg, read_1_end,
	      read_2_idx, read_2_beg, read_2_end);
	const BaseVec & bv1 = bvv[read_1_idx];
	const BaseVec & bv2 = bvv[read_2_idx];
	assert(read_1_beg < bv1.size());
	assert(read_1_end < bv1.size());
	assert(read_2_beg < bv2.size());
	assert(read_2_end < bv2.size());
	assert(read_1_beg != read_1_end);
	assert(read_2_beg != read_2_end);

	bool is_rc_1 = (read_1_beg > read_1_end);
	bool is_rc_2 = (read_2_beg > read_2_end);

	if (read_1_end < read_1_beg)
		std::swap(read_1_beg, read_1_end);
	if (read_2_end < read_2_beg)
		std::swap(read_2_beg, read_2_end);

	unsigned long len_1, len_2;
	len_1 = read_1_end - read_1_beg + 1;
	len_2 = read_2_end - read_2_beg + 1;
	assert(len_1 == len_2);
	assert(len_1 >= min_overlap_len);
	assert_seed_valid(bv1, bv2, read_1_beg, read_2_beg, len_1,
			  is_rc_1, is_rc_2, "OVERLAP");
}

static void extend_seed(const BaseVec & bv1,
			const BaseVec & bv2,
			unsigned & pos1,
			unsigned & pos2,
			unsigned & len,
			const bool is_rc1,
			const bool is_rc2)
{
	if (is_rc1 == is_rc2) {
		unsigned max_left_extend = std::min(pos1, pos2);
		unsigned left_extend = 0;
		while (left_extend < max_left_extend) {
			if (bv1[pos1 - (left_extend + 1)] == bv2[pos2 - (left_extend + 1)])
				left_extend++;
			else
				break;
		}
		unsigned max_right_extend = std::min(bv1.size() - (pos1 + len),
						     bv2.size() - (pos2 + len));
		unsigned right_extend = 0;
		while (right_extend < max_right_extend) {
			if (bv1[pos1 + len + right_extend] == bv2[pos2 + len + right_extend])
				right_extend++;
			else
				break;
		}
		len += left_extend + right_extend;
		pos1 -= left_extend;
		pos2 -= left_extend;
	} else {
		assert(!is_rc1 && is_rc2);

		unsigned max_left_extend = std::min(pos1, bv2.size() - (pos2 + len));
		unsigned left_extend = 0;
		while (left_extend < max_left_extend) {
			if (bv1[pos1 - (left_extend + 1)] ==
			    (3 ^ bv2[pos2 + len + left_extend]))
				left_extend++;
			else
				break;
		}
		unsigned max_right_extend = std::min(bv1.size() - (pos1 + len), pos2);
		unsigned right_extend = 0;
		while (right_extend < max_right_extend) {
			if (bv1[pos1 + len + right_extend] ==
			    (3 ^ bv2[pos2 - (right_extend + 1)]))
				right_extend++;
			else
				break;
		}
		len += left_extend + right_extend;
		pos1 -= left_extend;
		pos2 -= right_extend;
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
	const BaseVec & bv1 = bvv[occ1.get_read_id()];
	const BaseVec & bv2 = bvv[occ2.get_read_id()];
	unsigned pos1 = occ1.get_read_pos();
	unsigned pos2 = occ2.get_read_pos();
	unsigned len = K;
	const bool is_rc1 = occ1.is_rc();
	const bool is_rc2 = occ2.is_rc();
	assert_seed_valid(bv1, bv2, pos1, pos2, len, is_rc1, is_rc2);
	extend_seed(bv1, bv2, pos1, pos2, len, is_rc1, is_rc2);
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
static void
overlaps_from_kmer_seed(const std::vector<KmerOccurrence> & occs,
			const BaseVecVec &bvv,
			const unsigned min_overlap_len,
			const unsigned max_edits,
			OverlapVecVec &ovv,
			unsigned long & num_overlaps,
			unsigned long & num_pairs_considered)
{
	Overlap o;
	for (size_t i = 0; i < occs.size(); i++) {
		for (size_t j = i + 1; j < occs.size(); j++) {
			num_pairs_considered++;
			KmerOccurrence occ1 = occs[i];
			KmerOccurrence occ2 = occs[j];
			if (occ1.is_rc() && !occ2.is_rc())
				std::swap(occ1, occ2);
			if (!find_overlap(bvv, occ1, occ2,
				          min_overlap_len, max_edits, K, o))
				continue;

			OverlapVecVec::OverlapSet & os = ovv[occ1.get_read_id()];
			if (os.find(o) != os.end())
				continue;

			assert_overlap_valid(o, bvv,
					     min_overlap_len, max_edits);
			os.insert(o);
			num_overlaps++;
		}
	}
}

template <unsigned K>
static void load_kmer_occurrences(const BaseVecVec &bvv,
				  std::unordered_map<Kmer<K>,
				  		     std::vector<KmerOccurrence> > &occ_map)
{
	info("Finding all occurrences of %u-mers in the reads", K);
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

			for (size_t k = 0; k < K - 1; k++)
				assert(fwd_kmer[k] == bv[k + j - (K - 1)]);
			for (size_t k = 0; k < K - 1; k++)
				assert(rev_kmer[K - 1 - k] == (3 ^ bv[k + j - (K - 1)]));

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
	info("Loaded %lu %u-mer occurrences into hash map",
	     num_kmer_occurrences, K);
}

template <unsigned K>
static void compute_overlaps(const BaseVecVec &bvv, 
			     const unsigned min_overlap_len,
			     const unsigned max_edits,
			     OverlapVecVec &ovv)
{
	typedef std::unordered_map<Kmer<K>, std::vector<KmerOccurrence> >
		KmerOccurrenceMap;

	ovv.clear();
	ovv.resize(bvv.size());

	KmerOccurrenceMap occ_map;

	load_kmer_occurrences(bvv, occ_map);

	info("Finding overlaps from %u-mer seeds", K);
	unsigned long num_overlaps = 0;
	unsigned long num_pairs_considered = 0;
	for (auto kmer_occ_pair : occ_map) {
		overlaps_from_kmer_seed<K>(kmer_occ_pair.second, bvv,
					   min_overlap_len, max_edits,
					   ovv, num_overlaps,
					   num_pairs_considered);
	}
	info("Found %lu overlaps", num_overlaps);
	info("Considered %lu read pairs", num_pairs_considered);
}

static const char *optstring = "l:e:h";
static const struct option longopts[] = {
	{"min-overlap-len", required_argument, NULL, 'l'},
	{"max-edits",   required_argument, NULL, 'e'},
	END_LONGOPTS
};

DEFINE_USAGE(
"Usage: compute-overlaps BVV_FILE OVERLAPS_FILE\n"
"\n"
"  -l, --min-overlap-len=LEN\n"
"  -e, --max-edits=MAX_EDITS\n"
"  -h, --help\n"
);

int main(int argc, char *argv[])
{
	int c;
	unsigned min_overlap_len = 25;
	unsigned max_edits = 0;
	for_opt(c) {
		switch (c) {
		case 'l':
			min_overlap_len = parse_long(optarg, "--min-overlap-len",
						     16, UINT_MAX);
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

	info("Loading reads from \"%s\"", argv[0]);
	BaseVecVec bvv(argv[0]);
	info("Loaded %zu reads from \"%s\"", bvv.size(), argv[0]);
	OverlapVecVec ovv;

	if (min_overlap_len < 24)
		compute_overlaps<16>(bvv, min_overlap_len, max_edits, ovv);
	else if (min_overlap_len < 32)
		compute_overlaps<24>(bvv, min_overlap_len, max_edits, ovv);
	else if (min_overlap_len < 40)
		compute_overlaps<32>(bvv, min_overlap_len, max_edits, ovv);
	else if (min_overlap_len < 48)
		compute_overlaps<40>(bvv, min_overlap_len, max_edits, ovv);
	else
		compute_overlaps<48>(bvv, min_overlap_len, max_edits, ovv);

	info("Writing overlaps to \"%s\"", argv[1]);
	ovv.write(argv[1]);
	info("Done writing \"%s\"", argv[1]);
}
