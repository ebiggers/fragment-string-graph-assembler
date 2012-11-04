#include "Overlap.h"
#include "BaseVecVec.h"
#include "DirectedStringGraph.h"

DEFINE_USAGE(
"Usage: build-graph READS_FILE OVERLAPS_FILE GRAPH_FILE\n"
);

template <typename GRAPH_t>
static void add_edge_from_overlap(const BaseVecVec & bvv, const Overlap & o,
				  GRAPH_t & graph)
{
	unsigned long f_idx;
	unsigned long f_beg;
	unsigned long f_end;
	unsigned long g_idx;
	unsigned long g_beg;
	unsigned long g_end;

	o.get(f_idx, f_beg, f_end, g_idx, g_beg, g_end);

	const BaseVec & f = bvv[f_idx];
	const BaseVec & g = bvv[g_idx];

	// Must not be a contained overlap
	assert(!((f_beg == 0 && f_end == f.size() - 1)
	    || (f_beg == f.size() - 1 && f_end == 0)
	    || (g_beg == 0 && g_end == g.size() - 1)
	    || (g_beg == g.size() - 1 && g_end == 0)));

	if (f_beg > f_end) {
		std::swap(f_idx, g_idx);
		std::swap(f_beg, g_beg);
		std::swap(f_end, g_end);
	}
	if (f_beg > 0) {
		if (g_beg < g_end) {
			/*
			 *  f.B --------------> f.E
			 *         g.B ----------------> g.E
			 *
			 *  Add f.E => g.E, g.B => f.B
			 *
			 *  Or bidirected edge:
			 *  
			 *  f >----------> g
			 *    
			 *     f.E -> g.E label: g[g_end + 1 ... g.size() - 1]
			 *     g.B -> f.B label: f[0 ... f_beg - 1]
			 */

			graph.add_edge_pair(f_idx, g_idx,
					    graph.TAG_F_E | graph.TAG_G_E,
					    g, g_end + 1, g.size() - 1,
					    f, f_beg - 1, 0);
		} else {

			/*
			 *  f.B --------------> f.E
			 *         g.E <---------------  g.B
			 *
			 *  Add f.E => g.B, g.E => f.B
			 *
			 *  Or bidirected edge:
			 *  
			 *  f <----------> g
			 *    
			 *     f.E -> g.B label: g[g_end - 1 ... 0]
			 *     g.E -> f.B label: f[f_beg - 1 ... 0]
			 */

			graph.add_edge_pair(f_idx, g_idx,
					    graph.TAG_F_E | graph.TAG_G_B,
					    g, g_end - 1, 0,
					    f, f_beg - 1, 0);
		}
	} else {
		if (g_beg < g_end) {

			/*
			 *        f.B ---------------> f.E
			 * g.B --------------> g.E
			 *
			 *  Add f.B => g.B, g.E => f.E
			 *
			 *  Or bidirected edge:
			 *  
			 *  f <----------< g
			 *    
			 *     f.B -> g.B label: g[g_beg - 1 ... 0]
			 *     g.E -> f.E label: f[f_end + 1 ... f.size() - 1]
			 */

			graph.add_edge_pair(f_idx, g_idx,
					    graph.TAG_F_B | graph.TAG_G_B,
					    g, g_beg - 1, 0,
					    f, f_beg + 1, f.size() - 1);
		} else {

			/*
			 *        f.B ---------------> f.E
			 * g.E <-------------- g.B
			 *
			 *  Add f.B => g.E, g.B => f.E
			 *
			 *  Or bidirected edge:
			 *  
			 *  f <----------> g
			 *    
			 *     f.B -> g.E label: g[g_beg + 1 ... g.size() - 1]
			 *     g.B -> f.E label: f[f_end + 1 ... f.size() - 1]
			 */
			graph.add_edge_pair(f_idx, g_idx,
					    graph.TAG_F_B | graph.TAG_G_E,
					    g, g_beg + 1, g.size() - 1,
					    f, f_end + 1, f.size() - 1);
		}
	}
}

template <typename GRAPH_t>
static void build_graph(const BaseVecVec & bvv, const OverlapVecVec & ovv,
			GRAPH_t & graph)
{
	for (auto overlap_set : ovv) {
		for (const Overlap & o : overlap_set) {
			assert_overlap_valid(o, bvv, 0, 0);
			add_edge_from_overlap(bvv, o, graph);
		}
	}
	info("String graph has %zu vertices and %zu edges",
	     graph.num_vertices(), graph.num_edges());
	info("Average of %.2f edges per vertex",
	     graph.num_edges() ?
	     	double(graph.num_edges()) / double(graph.num_vertices()) : 0);
}

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 4);
	const char *reads_file = argv[1];
	const char *overlaps_file = argv[2];
	const char *graph_file = argv[3];

	info("Reading reads from \"%s\"", reads_file);
	BaseVecVec bvv(reads_file);
	info("Loaded %zu reads from \"%s\"", bvv.size(), reads_file);

	info("Loading overlaps from \"%s\"", overlaps_file);
	OverlapVecVec ovv(overlaps_file);
	info("Done loading overlaps");

	assert(ovv.size() == bvv.size());

	DirectedStringGraph graph(bvv.size());

	info("Building string graph from overlaps");
	build_graph(bvv, ovv, graph);

	info("Writing string graph to \"%s\"", graph_file);
	graph.write(graph_file);
	info("Done");
}
