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

	// Per-vertex mark.  Initially set to VACANT.
	std::vector<unsigned char> vertex_marks(graph.num_vertices(), VACANT);

	// Per-edge boolean indicating whether each edge is to be deleted as
	// part of the transitive reduction or not.
	std::vector<bool> reduce_edge(graph.num_edges(), false);

	// Iterate through every vertex @v in the graph that has outgoing edges.
	for (size_t i = 0; i < vertices.size(); i++) {
		GraphVertex & v = vertices[i];
		const std::vector<unsigned long> & edge_indices = v.edge_indices();
		if (edge_indices.size() == 0)
			continue;

		// Mark each vertex adjacent to @v as INPLAY.
		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			assert(e.get_v1_idx() == i);
			assert(e.get_v2_idx() != i);
			vertex_marks[e.get_v2_idx()] = INPLAY;
		}

		// Length of the longest sequence label on the edges leaving
		// vertex @v.
		size_t longest = edges[edge_indices.back()].get_seq().size();

		// For each outgoing edge from v => w in order of labeled
		// sequence length, consider each vertex w that is still marked
		// INPLAY.
		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			unsigned long w_idx = e.get_v2_idx();
			if (vertex_marks[w_idx] == INPLAY) {
				// e must be an irreducible edge if w is still
				// marked INPLAY at this point, since all
				// shorter edges were already considered.
				//
				// Now, consider the edges leaving vertex w.
				// Each such edge that goes to a vertex marked
				// INPLAY must be directly reachable from v, and
				// therefore the edge must be removed.
				GraphVertex & w = vertices[w_idx];
				const std::vector<unsigned long> & w_edge_indices = w.edge_indices();
				for (size_t k = 0; k < w_edge_indices.size(); k++) {
					GraphEdge & e2 = edges[w_edge_indices[k]];
					if (e2.get_seq().size() <= longest) {
						if (vertex_marks[w_idx] == INPLAY) {
							vertex_marks[w_idx] = ELIMINATED;
						}
					}
				}
			}
		}

		#if 0
		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			unsigned long w_idx = e.get_v2_idx();
			GraphVertex & w = vertices[w_idx];
			const std::vector<unsigned long> & w_edge_indices = w.edge_indices();
			for (size_t k = 0; k < w_edge_indices.size(); k++) {
				if (k == 0) { // TODO: fuzz parameter
					if (vertex_marks[w_idx] == INPLAY) {
						vertex_marks[w_idx] = ELIMINATED;
					}
				}
			}
		}
		#endif

		// Once again, go through the outgoing edges from v.  For each
		// vertex marked ELIMINATED, mark the corresponding edge for
		// reduction.  Return both INPLAY and ELIMINATED vertices to
		// VACANT status.
		for (size_t j = 0; j < edge_indices.size(); j++) {
			GraphEdge & e = edges[edge_indices[j]];
			unsigned long w_idx = e.get_v2_idx();
			if (vertex_marks[w_idx] == ELIMINATED) {
				reduce_edge[edge_indices[j]] = true;
			}
			vertex_marks[w_idx] = VACANT;
		}
	}

	info("Transitive reduction algorithm complete.  Now updating the strign graph");

	// Modify the graph.
	std::vector<unsigned long> new_edge_indices(edges.size());
	size_t i, j;
	for (i = 0, j = 0; i < edges.size(); i++) {
		if (!reduce_edge[i]) {
			new_edge_indices[i] = j;
			edges[j] = edges[i];
			j++;
		} else {
			new_edge_indices[i] = std::numeric_limits<unsigned long>::max();
		}
	}
	size_t num_original_edges = edges.size();
	size_t num_removed_edges = num_original_edges - j;
	edges.resize(j);

	info("Removing %zu of %zu edges (%.2f%%)",
	     num_removed_edges, num_original_edges,
	     edges.size() ? 100 * double(num_removed_edges) / num_original_edges : 0.0);
	for (GraphVertex & v : vertices) {
		std::vector<unsigned long> & edge_indices = v.edge_indices();
		for (i = 0, j = 0; i < edge_indices.size(); i++) {
			unsigned long new_edge_idx = new_edge_indices[edge_indices[i]];
			if (new_edge_idx != std::numeric_limits<unsigned long>::max()) {
				edge_indices[j++] = new_edge_idx;
			}
		}
		edge_indices.resize(j);
	}

	info("Done removing transitive edges");
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
