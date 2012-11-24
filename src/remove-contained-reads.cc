#include "Overlap.h"
#include "BaseVecVec.h"
#include "util.h"

DEFINE_USAGE(
"Usage: remove-contained-reads READS_FILE UNCONTAINED_READS_FILE\n"
"                              OVERLAPS_FILE UNCONTAINED_OVERLAPS_FILE\n"
"                              OLD_TO_NEW_INDICES_FILE\n"
"\n"
"Given a set of reads and all overlaps that were computed from them, find all\n"
"reads that are fully contained by another read and discard them, along with\n"
"the corresponding overlaps.\n"
"\n"
"If there are identical reads, only one of each is kept.\n"
"\n"
"Input:\n"
"      READS_FILE:     The set of reads from which the overlaps were found.\n"
"      OVERLAPS_FILE:  The set of overlaps, computed from the reads in\n"
"                       READS_FILE.\n"
"\n"
"Output:\n"
"      UNCONTAINED_READS_FILE:    The set of reads, with contained reads\n"
"                                 removed.\n"
"      UNCONTAINED_OVERLAPS_FILE: The set of overlaps, with overlaps with\n"
"                                 contained reads removed.\n"
"      OLD_TO_NEW_INDICES_FILE:   A map from the old read indices to the new\n"
"                                 read indices.\n"
);


int main(int argc, char **argv)
{
	USAGE_IF(argc != 6);
	const char *reads_file = argv[1];
	const char *uncontained_reads_file = argv[2];
	const char *overlaps_file = argv[3];
	const char *uncontained_overlaps_file = argv[4];
	const char *old_to_new_indices_file = argv[5];

	info("Loading reads from \"%s\"", reads_file);
	BaseVecVec bvv(reads_file);

	info("Loaded %zu reads", bvv.size());

	info("Loading overlaps from \"%s\"", overlaps_file);
	OverlapVecVec ovv(overlaps_file);

	assert(bvv.size() == ovv.size());

	std::vector<bool> read_contained(bvv.size(), false);

	info("Searching for overlaps indicating contained reads");
	foreach (const OverlapVecVec::OverlapSet & overlap_set, ovv) {
		foreach(const Overlap & o, overlap_set) {
			Overlap::read_idx_t f_idx;
			Overlap::read_pos_t f_beg;
			Overlap::read_pos_t f_end;
			Overlap::read_idx_t g_idx;
			Overlap::read_pos_t g_beg;
			Overlap::read_pos_t g_end;
			bool rc;

			assert_overlap_valid(o, bvv, 1, 0);

			o.get(f_idx, f_beg, f_end, g_idx, g_beg, g_end, rc);
			const BaseVec & f = bvv[f_idx];
			const BaseVec & g = bvv[g_idx];

			if ((f_beg == 0 && f_end == f.size() - 1))
				read_contained[f_idx] = true;
			else if (g_beg == 0 && g_end == g.size() - 1)
				read_contained[g_idx] = true;
		}
	}

	info("Computing new read indices");
	std::vector<size_t> old_to_new_indices(bvv.size());
	BaseVecVec uncontained_bvv;
	size_t i, j;
	for (i = 0, j = 0; i < bvv.size(); i++) {
		if (read_contained[i]) {
			old_to_new_indices[i] = std::numeric_limits<size_t>::max();
			bvv[i].destroy();
		} else {
			old_to_new_indices[i] = j++;
			uncontained_bvv.push_back(bvv[i]);
		}
	}
	size_t num_contained_reads = bvv.size() - uncontained_bvv.size();
	info("%zu of %zu reads were contained (%.2f%%)",
	     num_contained_reads, bvv.size(),
	     TO_PERCENT(num_contained_reads, bvv.size()));

	bvv.clear();

	info("Deleting overlaps for the contained reads");
	unsigned long num_overlaps_deleted = 0;
	unsigned long num_overlaps = 0;
	for (i = 0, j = 0; i < ovv.size(); i++) {
		if (!read_contained[i]) {
			std::set<Overlap> new_set;
			foreach(const Overlap & o, ovv[i]) {
				Overlap::read_idx_t f_idx;
				Overlap::read_idx_t g_idx;
				num_overlaps++;
				o.get_indices(f_idx, g_idx);
				if (!read_contained[f_idx] && !read_contained[g_idx]) {
					f_idx = old_to_new_indices[f_idx];
					g_idx = old_to_new_indices[g_idx];
					Overlap new_o(o);
					new_o.set_indices(f_idx, g_idx);
					new_set.insert(new_o);
				} else {
					num_overlaps_deleted++;
				}
			}
			ovv[j++] = new_set;
		}
	}
	info("Deleted %lu of %lu overlaps (%.2f%%)",
	     num_overlaps_deleted, num_overlaps,
	     TO_PERCENT(num_overlaps_deleted, num_overlaps));
	ovv.resize(j);

	assert(ovv.size() == uncontained_bvv.size());

	info("Writing uncontained reads to \"%s\"", uncontained_reads_file);
	uncontained_bvv.write(uncontained_reads_file);
	info("Writing uncontained overlaps to \"%s\"", uncontained_overlaps_file);
	ovv.write(uncontained_overlaps_file);

	info("Writing map from old read indices to new read indices to \"%s\"",
	     old_to_new_indices_file);
	{
		std::ofstream out(old_to_new_indices_file);
		boost::archive::binary_oarchive ar(out);
		ar << old_to_new_indices;
		out.close();
		if (!out)
			fatal_error_with_errno("Error writing to \"%s\"",
					       old_to_new_indices_file);
	}

	info("Done");
	return 0;
}
