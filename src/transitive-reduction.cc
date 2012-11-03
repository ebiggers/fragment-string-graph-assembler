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

	static const unsigned char VACANT = 0;
	static const unsigned char INPLAY = 1;
	static const unsigned char ELIMINATED = 2;

	std::vector<unsigned char> vertex_marks(graph.num_vertices(), VACANT);
	std::vector<bool> reduce_edge(graph.num_edges(), false);

	for (size_t i = 0; i < vertices.size(); i++) {
		GraphVertex & v = vertices[i];
		const std::vector<unsigned long> & edge_indices = v.edge_indices();
		if (edge_indices.size() == 0)
			continue;
		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			assert(e.get_v1_idx() == i);
			assert(e.get_v2_idx() != i);
			vertex_marks[e.get_v2_idx()] = INPLAY;
		}
		size_t longest = edges[edge_indices.back()].get_seq().size();
		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			unsigned long v2_idx = e.get_v2_idx();
			if (vertex_marks[v2_idx] == INPLAY) {
				GraphVertex & v2 = vertices[v2_idx];
				const std::vector<unsigned long> & v2_edge_indices = v2.edge_indices();
				for (size_t k = 0; k < v2_edge_indices.size(); k++) {
					GraphEdge & e = edges[v2_edge_indices[k]];
					if (e.get_seq().size() <= longest) {
						if (vertex_marks[v2_idx] == INPLAY) {
							vertex_marks[v2_idx] = ELIMINATED;
						}
					}
				}
			}
		}

		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			unsigned long v2_idx = e.get_v2_idx();
			GraphVertex & v2 = vertices[v2_idx];
			const std::vector<unsigned long> & v2_edge_indices = v2.edge_indices();
			for (size_t k = 0; k < v2_edge_indices.size(); k++) {
				if (k == 0) { // TODO: fuzz parameter
					if (vertex_marks[v2_idx] == INPLAY) {
						vertex_marks[v2_idx] = ELIMINATED;
					}
				}
			}
		}

		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			unsigned long v2_idx = e.get_v2_idx();
			if (vertex_marks[v2_idx] == ELIMINATED) {
				reduce_edge[edge_indices[j]] = true;
			}
			vertex_marks[v2_idx] = VACANT;
		}
	}
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
