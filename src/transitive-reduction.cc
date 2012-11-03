#include "Graph.h"
#include <algorithm>

DEFINE_USAGE(
"Usage: transitive-reduction GRAPH_FILE OUT_GRAPH_FILE"
);

class sort_by_edge_length {
private:
	const std::vector<GraphEdge> & _edges;
public:
	sort_by_edge_length(const std::vector<GraphEdge> & edges)
		: _edges(edges) { }
	bool operator()(unsigned long edge_idx_1, unsigned long edge_idx_2) {
		size_t edge_1_len = _edges[edge_idx_1].get_seq().size();
		size_t edge_2_len = _edges[edge_idx_2].get_seq().size();
		return (edge_1_len < edge_2_len);
	}
};

static void transitive_reduction(Graph & graph)
{
	info("Performing transitive reduction on string graph with "
	     "%zu vertices and %zu edges", graph.num_vertices(), graph.num_edges());

	std::vector<GraphVertex> & vertices = graph.vertices();
	std::vector<GraphEdge> & edges = graph.edges();

	sort_by_edge_length cmp(edges);
	info("Sorting adjacency lists of vertices by edge length");
	for (GraphVertex & v : vertices)
		std::sort(v.edge_indices().begin(), v.edge_indices().end(), cmp);

	std::vector<bool> vertex_vacant(graph.num_vertices(), true);
	std::vector<bool> reduce_edge(graph.num_edges(), false);

	info("Done performing transitive reduction");
}

int main(int argc, char **argv)
{
	USAGE_IF(argc != 3);
	info("Loading string graph from \"%s\"", argv[1]);
	Graph graph(argv[1]);
	transitive_reduction(graph);
	info("Writing string graph to \"%s\"", argv[2]);
	graph.write(argv[2]);
}
