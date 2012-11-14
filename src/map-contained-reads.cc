#include "Overlap.h"
#include "BaseVecVec.h"
#include "AnyStringGraph.h"
#include "util.h"

DEFINE_USAGE(
"Usage: map-contained-reads ORIG_READS_FILE ORIG_OVERLAPS_FILE\n"
"                           OLD_NEW_INDICES_MAP_FILE GRAPH_FILE OUT_GRAPH_FILE\n"
"\n"
"Map contained reads back into a graph.\n"
"\n"
"Input:\n"
"      READS_FILE:     The set of reads from which the overlaps were found.\n"
"      OVERLAPS_FILE:  The set of overlaps, computed from the reads in\n"
"                       READS_FILE.\n"
"      OLD_NEW_INDICES_MAP_FILE:  A map from the old read indices to the new\n"
"                                 read indices.\n"
"      GRAPH_FILE:\n"
"\n"
"Output:\n"
"      OUT_GRAPH_FILE:\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 6);
	const char * const orig_reads_file            = argv[1];
	const char * const orig_overlaps_file        = argv[2];
	const char * const old_new_indices_map_file  = argv[3];
	const char * const graph_file                = argv[4];
	const char * const out_graph_file            = argv[5];

	BaseVecVec orig_reads(orig_reads_file);
	OverlapVecVec orig_overlaps(orig_overlaps_file);

	std::vector<size_t> old_to_new_indices;
	{
		std::ifstream is(old_new_indices_map_file);
		boost::archive::binary_iarchive ar(is);
		ar >> old_to_new_indices;
	}


	// Build a map from original indices to contained read indices
	assert(old_to_new_indices.size() == orig_reads.size()
	       && orig_overlaps.size() == orig_reads.size());

	std::vector<size_t> contained_read_indices;
	std::vector<size_t> orig_to_contained_indices(orig_reads.size());
	std::vector<size_t> new_to_old_indices(old_to_new_indices.size(), ~size_t(0));
	for (size_t i = 0; i < orig_reads.size(); i++) {
		new_to_old_indices[old_to_new_indices[i]] = i;
		if (old_to_new_indices[i] == ~size_t(0)) {
			orig_to_contained_indices[i] = contained_read_indices.size();
			contained_read_indices.push_back(i);
		} else {
			orig_to_contained_indices[i] = ~size_t(0);
		}
	}

	size_t num_contained_reads = contained_read_indices.size();

	std::vector<std::vector<Overlap *>> containing_overlaps(num_contained_reads);

	info("Examining %zu contained reads (%.2f%% of total)",
	     num_contained_reads,
	     orig_reads.size() ?
	     	num_contained_reads * 100 / double(orig_reads.size()) : 0);

	// Find all overlaps involving each contained read
	size_t num_contained_overlaps = 0;
	foreach (const std::set<Overlap> & overlap_set, orig_overlaps) {
		foreach(const Overlap & o, overlap_set) {
			Overlap::read_idx_t f_idx;
			Overlap::read_pos_t f_beg;
			Overlap::read_pos_t f_end;
			Overlap::read_idx_t g_idx;
			Overlap::read_pos_t g_beg;
			Overlap::read_pos_t g_end;
			bool rc;

			assert_overlap_valid(o, orig_reads, 0, 0);

			o.get(f_idx, f_beg, f_end, g_idx, g_beg, g_end, rc);
			const BaseVec & f = orig_reads[f_idx];
			const BaseVec & g = orig_reads[g_idx];

			size_t idx;
			if ((f_beg == 0 && f_end == f.size() - 1)) {
				if (g_beg == 0 && g_end == g.size() - 1)
					idx = std::min(f_idx, g_idx);
				else
					idx = f_idx;

				assert(old_to_new_indices[idx] == ~size_t(0));
				idx = orig_to_contained_indices[idx];
				assert(idx != ~size_t(0));
				containing_overlaps[idx].push_back(&o);
				num_contained_overlaps++;
			} else if (g_beg == 0 && g_end == g.size() - 1) {
				idx = g_idx;
				assert(old_to_new_indices[idx] == ~size_t(0));
				idx = orig_to_contained_indices[idx];
				assert(idx != ~size_t(0));
				containing_overlaps[idx].push_back(&o);
				num_contained_overlaps++;
			}
		}
	}

	info("%zu overlaps are contained (%.2f per contained read on average)",
	     num_contained_overlaps, 
	     (num_contained_reads ?
	     	double(num_contained_overlaps) / num_contained_reads : 0));

	std::vector<Overlap *> shortest_overhang_overlaps(num_contained_reads, NULL);
	std::vector<Overlap::read_pos_t> overhang_lens(num_contained_reads,
						       std::numeric_limits<Overlap::read_pos_t>::max());

	// Find the shortest overhanging overlap for each contained read
	for (size_t i = 0; i < num_contained_reads; i++) {
		size_t orig_read_idx = contained_read_indices[i];

		assert(containing_overlaps[i].size() != 0);

		foreach (Overlap * o, containing_overlaps[i]) {
			Overlap::read_idx_t f_idx;
			Overlap::read_pos_t f_beg;
			Overlap::read_pos_t f_end;
			Overlap::read_idx_t g_idx;
			Overlap::read_pos_t g_beg;
			Overlap::read_pos_t g_end;
			bool rc;

			o->get(f_idx, f_beg, f_end, g_idx, g_beg, g_end, rc);


			if (orig_read_idx == g_idx) {
				std::swap(f_idx, g_idx);
				std::swap(f_beg, g_beg);
				std::swap(f_end, g_end);
			}

			const BaseVec & f = orig_reads[f_idx];
			const BaseVec & g = orig_reads[g_idx];

			Overlap::read_pos_t overhang_len;
			assert(g.length() >= g_end + 1);
			overhang_len = g.length() - 1 - g_end;

			if (overhang_len < overhang_lens[i]) {
				shortest_overhang_overlaps[i] = o;
				overhang_lens[i] = overhang_len;
			}
		}
		assert(shortest_overhang_overlaps[i] != NULL);
	}


	AnyStringGraph graph(graph_file);

	// Map the contained reads into the graph, given the overlap to use to
	// map each contained read.
	for (size_t i = 0; i < num_contained_reads; i++) {
		size_t orig_read_idx = contained_read_indices[i];
		Overlap * o = shortest_overhang_overlaps[i];
		Overlap::read_idx_t f_idx;
		Overlap::read_idx_t g_idx;
		o->get_indices(f_idx, g_idx);
		assert(new_to_old_indices[orig_read_idx] != ~size_t(0));
		size_t new_read_idx = new_to_old_indices[orig_read_idx];
		assert(new_to_old_indices[f_idx] != ~size_t(0));
		assert(new_to_old_indices[g_idx] != ~size_t(0));
		f_idx = new_to_old_indices[f_idx];
		g_idx = new_to_old_indices[g_idx];
		o->set_indices(f_idx, g_idx);
		graph.map_contained_read(new_read_idx, *o);
	}

	return 0;
}
