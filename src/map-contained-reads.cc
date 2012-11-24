#include "Overlap.h"
#include "BaseVecVec.h"
#include "AnyStringGraph.h"
#include "util.h"

DEFINE_USAGE(
"Usage: map-contained-reads ORIG_READS_FILE ORIG_OVERLAPS_FILE\n"
"                           OLD_TO_NEW_INDICES_FILE GRAPH_FILE OUT_GRAPH_FILE\n"
"\n"
"Map contained reads back into a graph.\n"
"\n"
"Input:\n"
"      READS_FILE:     The set of reads from which the overlaps were found.\n"
"      OVERLAPS_FILE:  The set of overlaps, computed from the reads in\n"
"                       READS_FILE.\n"
"      OLD_TO_NEW_INDICES_FILE:   A map from the old read indices to the new\n"
"                                 read indices.\n"
"      GRAPH_FILE:     The graph to map the contained reads into.\n"
"\n"
"Output:\n"
"      OUT_GRAPH_FILE: The output graph into which the contained reads have.\n"
"                      been mapped.\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 6);
	const char * const orig_reads_file         = argv[1];
	const char * const orig_overlaps_file      = argv[2];
	const char * const old_to_new_indices_file = argv[3];
	const char * const graph_file              = argv[4];
	const char * const out_graph_file          = argv[5];

	info("Loading original reads from \"%s\"", orig_reads_file);
	BaseVecVec orig_reads(orig_reads_file);

	info("Loading original overlaps from \"%s\"", orig_overlaps_file);
	OverlapVecVec orig_overlaps(orig_overlaps_file);

	info("Loading map from old to new read indices from \"%s\"", old_to_new_indices_file);
	std::vector<size_t> old_to_new_indices;
	{
		std::ifstream is(old_to_new_indices_file);
		boost::archive::binary_iarchive ar(is);
		ar >> old_to_new_indices;
	}

	size_t num_orig_reads = orig_reads.size();

	assert(old_to_new_indices.size() == num_orig_reads);
	assert(orig_overlaps.size() == num_orig_reads);

	// Count the number of contained reads.
	size_t num_contained_reads = 0;
	for (size_t i = 0; i < num_orig_reads; i++)
		if (old_to_new_indices[i] == ~size_t(0))
			num_contained_reads++;

	// The number of uncontained reads is the number of original reads minus
	// the number of contained reads.
	size_t num_uncontained_reads = num_orig_reads - num_contained_reads;

	info("%zu of %zu original reads were contained (%.2f%%)",
	     num_contained_reads, num_orig_reads,
	     TO_PERCENT(num_contained_reads, num_orig_reads));

	// Use the old to new indices map to build the reverse map from the new
	// to old read indices.  At the same time, construct a sequential
	// numbering of the contained reads, and build a map from the original
	// read indices to contained read indices and a reverse map from the
	// contained read indices to the original read indices.
	std::vector<size_t> contained_to_old_indices(num_contained_reads, ~size_t(0));
	std::vector<size_t> new_to_old_indices(num_uncontained_reads, ~size_t(0));
	std::vector<size_t> old_to_contained_indices(num_orig_reads, ~size_t(0));
	size_t j = 0;
	for (size_t i = 0; i < orig_reads.size(); i++) {
		if (old_to_new_indices[i] == ~size_t(0)) {
			// The read exists in the original reads but not in the
			// new reads, so it is a contained read.
			old_to_contained_indices[i] = j;
			contained_to_old_indices[j++] = i;
		} else {
			// The read exists in both the original and new reads,
			// so it is not a contained read.  It may have been
			// renumbered.
			assert(old_to_new_indices[i] < num_uncontained_reads);
			assert(new_to_old_indices[old_to_new_indices[i]] == ~size_t(0));
			new_to_old_indices[old_to_new_indices[i]] = i;
		}
	}
	assert(j == num_contained_reads);

	// Find the shortest overhanging overlap for each contained read
	info("Finding the shortest overhanging overlap for each contained read");

	std::vector<const Overlap *>
		shortest_overhang_overlaps(num_contained_reads, NULL);

	std::vector<Overlap::read_pos_t>
		shortest_overhang_lens(num_contained_reads,
				       std::numeric_limits<Overlap::read_pos_t>::max());

	foreach (std::set<Overlap> & overlap_set, orig_overlaps) {
		foreach(const Overlap & o, overlap_set) {
			Overlap::read_idx_t f_idx;
			Overlap::read_pos_t f_beg;
			Overlap::read_pos_t f_end;
			Overlap::read_idx_t g_idx;
			Overlap::read_pos_t g_beg;
			Overlap::read_pos_t g_end;
			bool rc;

			Overlap::read_idx_t uncontained_read_orig_idx;
			Overlap::read_pos_t uncontained_read_overlap_beg;
			Overlap::read_pos_t uncontained_read_overlap_end;

			assert_overlap_valid(o, orig_reads, 1, 0);

			o.get(f_idx, f_beg, f_end, g_idx, g_beg, g_end, rc);
			const BaseVec & f = orig_reads[f_idx];
			const BaseVec & g = orig_reads[g_idx];
			const BaseVec * uncontained_read;

			assert(f_idx < num_orig_reads);
			assert(g_idx < num_orig_reads);

			size_t contained_read_orig_idx = ~size_t(0);
			if (f_beg == 0 && f_end == f.size() - 1) {

				// Read f is contained
				assert(old_to_new_indices[f_idx] == ~size_t(0));
				assert(old_to_contained_indices[f_idx] != ~size_t(0));

				// ... but only count this overlap if g is NOT
				// contained
				if (old_to_new_indices[g_idx] != ~size_t(0)) {
					contained_read_orig_idx      = f_idx;
					uncontained_read             = &g;
					uncontained_read_orig_idx    = g_idx;
					uncontained_read_overlap_beg = g_beg;
					uncontained_read_overlap_end = g_end;
				}
			} else if (g_beg == 0 && g_end == g.size() - 1) {
				// Read g is contained
				assert(old_to_new_indices[g_idx] == ~size_t(0));
				assert(old_to_contained_indices[g_idx] != ~size_t(0));

				// ... but only count this overlap is f is NOT
				// contained
				if (old_to_new_indices[f_idx] != ~size_t(0)) {
					contained_read_orig_idx      = g_idx;
					uncontained_read             = &f;
					uncontained_read_orig_idx    = f_idx;
					uncontained_read_overlap_beg = f_beg;
					uncontained_read_overlap_end = f_end;
				}
			}

			if (contained_read_orig_idx != ~size_t(0)) {
				Overlap::read_pos_t overhang_len;
				// This is a containing overlap, and the read
				// with original index @contained_read_orig_idx
				// is contained.  Compute the overhang length.
				if (rc) {
					//
					//
					//      ---------->
					//  <------------------------
					//
					//                 |overhang|
					overhang_len = uncontained_read_overlap_beg;
				} else {
					//
					//
					//      ---------->
					//  ------------------------>
					//
					//                 |overhang|
					assert(uncontained_read->length() >=
					       uncontained_read_overlap_end + 1);
					overhang_len = uncontained_read->length() -
						       (uncontained_read_overlap_end + 1);
				}

				// Index of this contained read in the
				// sequential numbering of the contained reads
				size_t contained_idx =
					old_to_contained_indices[contained_read_orig_idx];

				// Update the shortest overhang if this overhang
				// is shorter than the previous shortest one.
				if (overhang_len < shortest_overhang_lens[contained_idx]) {
					shortest_overhang_overlaps[contained_idx] = &o;
					shortest_overhang_lens[contained_idx] = overhang_len;
				}
			}
		}
	}

	info("Reading string graph from \"%s\"", graph_file);
	AnyStringGraph graph(graph_file);

	// Map the contained reads into the graph, given the overlap to use to
	// map each contained read.
	for (size_t i = 0; i < num_contained_reads; i++) {

		// Original index of this contained read
		size_t contained_read_orig_idx = contained_to_old_indices[i];

		// The overlap in which the overhang from the end of this
		// contained read to its containing read is the shortest
		const Overlap *o = shortest_overhang_overlaps[i];

		// Original read index must be valid, and the read must be
		// contained, and it must have at least 1 overlap with an
		// uncontained read.
		assert(contained_read_orig_idx < num_orig_reads);
		assert(old_to_new_indices[contained_read_orig_idx] == ~size_t(0));
		assert(o != NULL);

		Overlap::read_idx_t f_idx, g_idx;
		Overlap::read_idx_t uncontained_read_dir;
		Overlap::read_idx_t uncontained_read_new_idx;

		o->get_indices(f_idx, g_idx);

		if (f_idx == contained_read_orig_idx)
			uncontained_read_new_idx = old_to_new_indices[g_idx];
		else
			uncontained_read_new_idx = old_to_new_indices[f_idx];

		assert(uncontained_read_new_idx < num_uncontained_reads);

		uncontained_read_dir = (o->is_rc() ? 1 : 0);

		info("Mapping read %zu of %zu (rc = %d)",
		     i + 1, num_contained_reads, o->is_rc());

		graph.map_contained_read(uncontained_read_new_idx,
					 uncontained_read_dir,
					 shortest_overhang_lens[i]);
	}

	info("Writing string graph to \"%s\"", out_graph_file);
	graph.write(out_graph_file);
	return 0;
}
