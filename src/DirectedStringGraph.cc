#include "DirectedStringGraph.h"
#include <fstream>
#include <ostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

class cmp_by_edge_length {
private:
	const std::vector<DirectedStringGraphEdge> & _edges;
public:
	cmp_by_edge_length(const std::vector<DirectedStringGraphEdge> & edges)
		: _edges(edges) { }

	bool operator()(unsigned long edge_idx_1, unsigned long edge_idx_2) const
	{
		return _edges[edge_idx_1].length() < _edges[edge_idx_2].length();
	}
};

void DirectedStringGraph::transitive_reduction()
{
	info("Performing transitive reduction on directed string graph with "
	     "%zu vertices and %zu edges", this->num_vertices(), this->num_edges());

	std::vector<DirectedStringGraphVertex> & vertices = this->vertices();
	std::vector<DirectedStringGraphEdge> & edges = this->edges();

	cmp_by_edge_length cmp(edges);
	info("Sorting adjacency lists of vertices by edge length");
	for (DirectedStringGraphVertex & v : vertices)
		std::sort(v.edge_indices().begin(), v.edge_indices().end(), cmp);

	static const unsigned char VACANT = 0;
	static const unsigned char INPLAY = 1;
	static const unsigned char ELIMINATED = 2;

	// Per-vertex mark.  Initially set to VACANT.
	std::vector<unsigned char> vertex_marks(this->num_vertices(), VACANT);

	// Per-edge boolean indicating whether each edge is to be deleted as
	// part of the transitive reduction or not.
	std::vector<bool> reduce_edge(this->num_edges(), false);

	// Iterate through every vertex @v in the graph that has outgoing edges.
	for (size_t v_idx = 0; v_idx < vertices.size(); v_idx++) {
		const DirectedStringGraphVertex & v = vertices[v_idx];

		if (v.out_degree() == 0)
			continue;

		// Mark each vertex adjacent to @v as INPLAY.
		for (const unsigned long edge_idx : v.edge_indices()) {
			const DirectedStringGraphEdge & e = edges[edge_idx];
			assert(e.get_v1_idx() == v_idx);
			assert(e.get_v2_idx() != v_idx);
			vertex_marks[e.get_v2_idx()] = INPLAY;
		}

		// Length of the longest sequence label on the edges leaving
		// vertex @v.
		const size_t longest = edges[v.edge_indices().back()].length();

		// For each outgoing edge from v => w in order of labeled
		// sequence length, consider each vertex w that is still marked
		// INPLAY.
		for (const unsigned long edge_idx : v.edge_indices()) {
			const unsigned long w_idx = edges[edge_idx].get_v2_idx();
			if (vertex_marks[w_idx] == INPLAY) {
				// The edge v -> w must be an irreducible edge
				// if w is still marked INPLAY at this point,
				// since all shorter edges were already
				// considered.
				//
				// Now, consider the edges leaving vertex w.
				// Each such edge that goes to a vertex marked
				// INPLAY must be directly reachable from v, and
				// therefore the edge must be removed.
				const DirectedStringGraphVertex & w = vertices[w_idx];
				for (const unsigned long w_edge_idx : w.edge_indices()) {
					const DirectedStringGraphEdge & e2 = edges[w_edge_idx];
					if (e2.length() > longest)
						break;
					if (vertex_marks[e2.get_v2_idx()] == INPLAY) {
						vertex_marks[e2.get_v2_idx()] = ELIMINATED;
					}
				}
			}
		}

		#if 0
		for (size_t j = 0; j < edge_indices.size(); j++) {
			DirectedStringGraphEdge & e = edges[edge_indices[j]];
			unsigned long w_idx = e.get_v2_idx();
			DirectedStringGraphVertex & w = vertices[w_idx];
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
		// neighboring vertex marked ELIMINATED, mark the corresponding
		// edge(s) for reduction.  Return both INPLAY and ELIMINATED
		// vertices to VACANT status.
		for (const unsigned long edge_idx : v.edge_indices())
			if (vertex_marks[edges[edge_idx].get_v2_idx()] == ELIMINATED)
				reduce_edge[edge_idx] = true;

		for (const unsigned long edge_idx : v.edge_indices())
			vertex_marks[edges[edge_idx].get_v2_idx()] = VACANT;
	}

	info("Transitive reduction algorithm complete.  Now updating the string graph");

	// Modify the graph
	std::vector<unsigned long> new_edge_indices(edges.size());
	size_t i, j;
	for (i = 0, j = 0; i < edges.size(); i++) {
		if (reduce_edge[i]) {
			new_edge_indices[i] = std::numeric_limits<unsigned long>::max();
		} else {
			new_edge_indices[i] = j;
			edges[j++] = edges[i];
		}
	}
	const size_t num_original_edges = edges.size();
	const size_t num_removed_edges = num_original_edges - j;
	edges.resize(j);

	info("Removing %zu of %zu edges (%.2f%%)",
	     num_removed_edges, num_original_edges,
	     (num_original_edges ?
		100 * double(num_removed_edges) / num_original_edges : 0.0));

	for (DirectedStringGraphVertex & v : vertices) {
		std::vector<unsigned long> & edge_indices = v.edge_indices();
		size_t j = 0;
		for (const unsigned long edge_idx : edge_indices) {
			const unsigned long new_edge_idx = new_edge_indices[edge_idx];
			if (new_edge_idx != std::numeric_limits<unsigned long>::max()) {
				edge_indices[j++] = new_edge_idx;
			}
		}
		edge_indices.resize(j);
	}

	info("Done removing transitive edges");
}
