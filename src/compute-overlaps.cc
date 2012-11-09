#include "BaseVecVec.h"
#include "Overlap.h"
#include <getopt.h>
#include "Kmer.h"

#include <unordered_map>
#include <ostream>

// Stores the location of a k-mer in the read set.
class KmerOccurrence {
private:
	// Index of the read containing the k-mer sequence
	unsigned long _read_id	: 32;

	// Position of the k-mer within the read (0-indexed)
	unsigned long _read_pos	: 31;

	// Whether the canonical k-mer is reverse-complement
	unsigned long _rc	: 1;
public:

	KmerOccurrence(unsigned long read_id, unsigned long read_pos, bool rc)
		: _read_id(read_id), _read_pos(read_pos), _rc(rc)
	{ }

	unsigned long get_read_id() const { return _read_id; }

	void swap_reads(KmerOccurrence & other) {
		unsigned long tmp;
		tmp = _read_id;
		_read_id = other._read_id;
		other._read_id = tmp;
		tmp = _read_pos;
		_read_pos = other._read_pos;
		other._read_pos = tmp;
	}

	unsigned long get_read_pos() const { return _read_pos; }

	bool is_rc() const { return _rc; }
	void flip_rc() { _rc = !_rc; }

	friend std::ostream & operator<<(std::ostream & os, const KmerOccurrence & occ)
	{
		return os << "KmerOccurrence { _read_id: " << occ._read_id <<
			", _read_pos: " << occ._read_pos <<
			", _rc: " <<occ._rc << "}";
	}
};



//
// Given a seed (an exactly matching sequence of bases of length @len, allowing
// for either forward or reverse-complement sequence) in the reads @bv1 and
// @bv2, extend it as far as possible on the left and right.
//
static void extend_seed(const BaseVec & bv1,
			const BaseVec & bv2,
			unsigned & pos1,
			unsigned & pos2,
			unsigned & len,
			const bool is_rc_1,
			const bool is_rc_2)
{
	assert_seed_valid(bv1, bv2, pos1, pos2, len, is_rc_1, is_rc_2);
	if (is_rc_1 == is_rc_2) {
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
		assert(!is_rc_1 && is_rc_2);

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

//
// Looks for an overlap seeded at the k-mers at @occ1 and @occ2.
//
// Returns %true and fills in the Overlap @o if a valid overlap is found.
//
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
	Overlap::read_pos_t pos1 = occ1.get_read_pos();
	Overlap::read_pos_t pos2 = occ2.get_read_pos();
	Overlap::read_pos_t len = K;
	const bool is_rc_1 = occ1.is_rc();
	const bool is_rc_2 = occ2.is_rc();

	// Note: in this function we should NOT be given a pair of k-mer
	// occurrences that are (reverse-complement, forward).  They always must
	// be switched to (forward, reverse-complement) first.
	assert(!(is_rc_1 && !is_rc_2));

	extend_seed(bv1, bv2, pos1, pos2, len, is_rc_1, is_rc_2);

	// The overlap must be at least @min_overlap_len base pairs long.
	if (len < min_overlap_len)
		return false;

	Overlap::read_pos_t read_1_beg = pos1;
	Overlap::read_pos_t read_1_end = pos1 + len - 1;
	Overlap::read_pos_t read_2_beg = pos2;
	Overlap::read_pos_t read_2_end = pos2 + len - 1;

	unsigned num_extremes = 0;
	if (read_1_beg == 0)
		num_extremes++;
	if (read_1_end == bv1.size() - 1)
		num_extremes++;
	if (read_2_beg == 0)
		num_extremes++;
	if (read_2_end == bv2.size() - 1)
		num_extremes++;

	// The overlap must have at least 2 extreme points.
	if (num_extremes < 2)
		return false;

	// If the overlap is of a read with itself, the overlap must be between
	// distinct regions.
	if (&bv1 == &bv2 &&
	    read_1_beg == read_2_beg &&
	    read_1_end == read_2_end)
		return false;

	// Eliminate overlaps that are definitely false because of how the reads
	// are laid out.  For example, if two reads happen to share the same
	// prefix longer than the minimum overlap length, but the remainder of
	// both reads are different, this is a false overlap.
	//
	// If two reads share the same prefix, the only way the overlap can be
	// true is if there are 3 extreme points (or 4, which would be an
	// overlap of identical, possibly reverse-complement reads, which is
	// allowed at this point).
	//
	// The same applies to reads that share the same suffix.

	{
		Overlap::read_pos_t _read_2_beg, _read_2_end;
		if (!is_rc_1 && is_rc_2) {
			_read_2_beg = (bv2.size() - 1) - read_2_end;
			_read_2_end = (bv2.size() - 1) - read_2_beg;
		} else {
			_read_2_beg = read_2_beg;
			_read_2_end = read_2_end;
		}

		// Shared prefix?
		if (read_1_beg == 0 && _read_2_beg == 0 && num_extremes < 3)
			return false;

		// Shared suffix?
		if (read_1_end == bv1.size() - 1 &&
		    _read_2_end == bv2.size() - 1 && num_extremes < 3)
			return false;
	}

	o.set(occ1.get_read_id(), read_1_beg, read_1_end,
	      occ2.get_read_id(), read_2_beg, read_2_end,
	      (!is_rc_1 && is_rc_2));
	return true;
}

//
// Finds all the overlaps that can be seeded at the occurrences of the canonical
// k-mer that has occurrences in the vector @occs.  Non-duplicate overlaps are
// added to the vector @ovv.
//
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
	// Consider each pair of k-mer occurrences only one time
	// (i.e. start j at i + 1, not 0)
	for (size_t i = 0; i < occs.size(); i++) {
		for (size_t j = i + 1; j < occs.size(); j++) {
			num_pairs_considered++;
			KmerOccurrence occ1 = occs[i];
			KmerOccurrence occ2 = occs[j];

			// The first occurrence is always set to the one with
			// lower read ID.
			if (occ1.get_read_id() > occ2.get_read_id())
				occ1.swap_reads(occ2);

			// If one occurrence is reverse-complement and the other
			// is not, always consider the first occurrence to the
			// forward and the second occurrence to be
			// reverse-complement.
			if (occ1.is_rc() && !occ2.is_rc()) {
				occ1.flip_rc();
				occ2.flip_rc();
			}

			// It may be the case that occ1.get_read_id() ==
			// occ2.get_read_id().  This happens if the same k-mer
			// occurs multiple times in one read, and this case is
			// allowed, but the overlap will be discarded unless two
			// *different* parts of the read overlap each other.

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

//
// Fills in the hash table @occ_map with a list of occurrences of each k-mer
// that appears in the reads @bvv.
//
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

//
// Compute overlaps.
//
// @bvv:
// 	Vector of reads.
//
// @min_overlap_len:
// 	Minimum length for each overlap.
//
// @max_edits:
// 	(Unimplemented)
//
// @ovv:
// 	Vector, indexed by read-id, into which a set of Overlaps for each read
// 	will be stored.
//
// Templatized by K, the length of the k-mer seed used to find overlaps.
template <unsigned K>
static void compute_overlaps(const BaseVecVec &bvv,
			     const unsigned min_overlap_len,
			     const unsigned max_edits,
			     OverlapVecVec &ovv)
{
	typedef std::unordered_map<Kmer<K>, std::vector<KmerOccurrence> >
		KmerOccurrenceMap;

	if (max_edits > 0)
		unimplemented();

	if (bvv.size() > Overlap::MAX_READ_IDX + 1) {
		fatal_error("'class Overlap' only supports up to %zu reads",
			    Overlap::MAX_READ_IDX + 1);
	}
	for (const BaseVec & bv : bvv) {
		if (bv.size() > Overlap::MAX_READ_LEN + 1) {
			fatal_error("'class Overlap' only supports reads up to "
				    "%zu bp long", Overlap::MAX_READ_LEN + 1);
		}
	}

	ovv.clear();
	ovv.resize(bvv.size());

	KmerOccurrenceMap occ_map;

	load_kmer_occurrences(bvv, occ_map);

	info("Finding overlaps from %u-mer seeds", K);
	unsigned long num_overlaps = 0;
	unsigned long num_pairs_considered = 0;
	for (auto kmer_occs_pair : occ_map) {
		overlaps_from_kmer_seed<K>(kmer_occs_pair.second, bvv,
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
"Usage: compute-overlaps READS_FILE OVERLAPS_FILE\n"
"\n"
"Computes all overlaps between reads in a set of reads.\n"
"\n"
"Input:\n"
"     READS_FILE:  FASTQ, FASTA, or binary reads (BaseVecVec) file\n"
"                  containing the read set.\n"
"\n"
"Output:\n"
"     OVERLAPS_FILE:  File to write the overlaps to.\n"
"\n"
"Options:\n"
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
	else if (min_overlap_len < 64)
		compute_overlaps<48>(bvv, min_overlap_len, max_edits, ovv);
	else if (min_overlap_len < 96)
		compute_overlaps<64>(bvv, min_overlap_len, max_edits, ovv);
	else if (min_overlap_len < 128)
		compute_overlaps<96>(bvv, min_overlap_len, max_edits, ovv);
	else
		compute_overlaps<128>(bvv, min_overlap_len, max_edits, ovv);

	info("Writing overlaps to \"%s\"", argv[1]);
	ovv.write(argv[1]);
	info("Done writing \"%s\"", argv[1]);
}
