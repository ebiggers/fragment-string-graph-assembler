#include "Overlap.h"
#include "BaseVec.h"
#include "BaseVecVec.h"

//
// Assert that the bases of the read @bv1 beginning at index @pos1 exactly match
// the bases of the read @bv2 beginning at index @pos2, for @len bases, where
// the sequence of read 1 is reverse-complemented iff @is_rc1 is %true and the
// sequence of read 2 is reverse-complemented iff @is_rc2 is %true.
//
void assert_seed_valid(const BaseVec & bv1,
		       const BaseVec & bv2,
		       const unsigned pos1,
		       const unsigned pos2,
		       const unsigned len,
		       const bool is_rc1, const bool is_rc2,
		       const char *description)
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

//
// Checks to make sure an overlap was correctly computed.
//
void assert_overlap_valid(const Overlap & o, const BaseVecVec & bvv,
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
