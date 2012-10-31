#include "Overlap.h"
#include "BaseVecVec.h"
#include "util.h"

DEFINE_USAGE(
"Usage: build-graph READS_FILE OVERLAPS_FILE Graph_FILE\n"
);

class GraphEdge {
private:
	unsigned long vertex_1_id : 32;
	unsigned long vertex_2_id : 32;
	BaseVec sequence;
public:
};

class GraphVertex {
private:
	std::vector<GraphEdge> edges;
public:
};

class Graph {
private:
	std::vector<GraphVertex> vertices;
public:
	enum {
		READ_BEGIN,
		READ_END
	};

	void write(const char *filename) const
	{
		unimplemented();
	}

	void add_edge(const unsigned long start_read_id, const int start_read_end,
		      const unsigned long end_read_id, const int end_read_end,
		      const BaseVec & bv,
		      const unsigned long beg, const unsigned long end)
	{
		unimplemented();
	}
};

typedef std::vector<GraphVertex> graph;

static void add_edge_from_overlap(const BaseVecVec & bvv, const Overlap & o,
				  Graph & graph)
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

	if (f_beg > f_end) {
		std::swap(f_idx, g_idx);
		std::swap(f_beg, g_beg);
		std::swap(f_end, g_end);
	}
	std::cout << o << std::endl;
	if (f_beg > 0) {
		if (g_beg < g_end) {

			graph.add_edge(g_idx, Graph::READ_BEGIN,
				       f_idx, Graph::READ_BEGIN,
				       f, f_beg, 0);

			graph.add_edge(f_idx, Graph::READ_END,
			               g_idx, Graph::READ_END,
				       g, g_end, g.size() - 1);
		} else {

			graph.add_edge(g_idx, Graph::READ_END,
				       f_idx, Graph::READ_BEGIN,
				       f, f_beg, 0);

			graph.add_edge(f_idx, Graph::READ_END,
			               g_idx, Graph::READ_BEGIN,
				       g, g_end, 0);
		}
	} else {
		if (g_beg < g_end) {

			graph.add_edge(f_idx, Graph::READ_BEGIN,
				       g_idx, Graph::READ_BEGIN,
				       g, g_beg, 0);

			graph.add_edge(g_idx, Graph::READ_END,
			               f_idx, Graph::READ_END,
				       f, f_end, f.size() - 1);

		} else {

			graph.add_edge(f_idx, Graph::READ_BEGIN,
				       g_idx, Graph::READ_END,
				       g, g.size() - 1, g_beg);

			graph.add_edge(g_idx, Graph::READ_BEGIN,
			               f_idx, Graph::READ_END,
				       f, f_end, f.size() - 1);

		}
	}
}

static void build_graph(const BaseVecVec & bvv, const OverlapVecVec & ovv,
			Graph & graph)
{
	for (auto overlap_set : ovv) {
		for (const Overlap & o : overlap_set) {
			add_edge_from_overlap(bvv, o, graph);
		}
	}
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

	Graph graph;

	info("Building graph");
	build_graph(bvv, ovv, graph);

	info("Writing graph to \"%s\"", graph_file);
	graph.write(graph_file);
	info("Done");
}
