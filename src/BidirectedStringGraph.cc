#include "BidirectedStringGraph.h"

void BidirectedStringGraph::transitive_reduction()
{
	info("Performing transitive reduction on bidirected string graph with "
	     "%zu vertices and %zu edges", this->num_vertices(), this->num_edges());

	std::vector<BidirectedStringGraphVertex> & vertices = this->vertices();
	std::vector<BidirectedStringGraphEdge> & edges = this->edges();

	info("Sorting adjacency lists of vertices by edge length");
	this->sort_adjlists_by_edge_len();

	static const unsigned char VACANT         = 0x0;
	static const unsigned char INPLAY_INWARD  = 0x1;
	static const unsigned char INPLAY_OUTWARD = 0x2;
	static const unsigned char ELIMINATED     = 0x4;
	static const unsigned char INPLAY         = INPLAY_INWARD | INPLAY_OUTWARD;

	// Per-vertex mark.  Initially set to VACANT.
	std::vector<unsigned char> vertex_marks(this->num_vertices(), VACANT);

	// Per-edge boolean indicating whether each edge is to be deleted as
	// part of the transitive reduction or not.
	std::vector<bool> reduce_edge(this->num_edges(), false);

	for (size_t v_idx = 0; v_idx < vertices.size(); v_idx++) {
		bool v_head_outward = true;
		do {
			const BidirectedStringGraphVertex & v = vertices[v_idx];

			BaseVec::size_type longest = 0;

			for (const edge_idx_t edge_idx : v.edge_indices()) {
				const BidirectedStringGraphEdge & e = edges[edge_idx];
				if (e.this_v_outward(v_idx) == v_head_outward) {
					if (e.other_v_outward(v_idx))
						vertex_marks[e.get_other_v_idx(v_idx)] = INPLAY_OUTWARD;
					else
						vertex_marks[e.get_other_v_idx(v_idx)] = INPLAY_INWARD;
					longest = std::max(longest, e.length());
				}
			}

			if (longest == 0)
				goto cont;

			for (const edge_idx_t edge_idx : v.edge_indices()) {
				const edge_idx_t w_idx = edges[edge_idx].get_other_v_idx(v_idx);
				if (vertex_marks[w_idx] & INPLAY) {
					const bool w_tail_outward =
							edges[edge_idx].other_v_outward(v_idx);
					const BidirectedStringGraphVertex & w = vertices[w_idx];
					for (const edge_idx_t w_edge_idx : w.edge_indices()) {
						const BidirectedStringGraphEdge & e2 = edges[w_edge_idx];
						if (e2.length() > longest)
							break;
						if (e2.this_v_outward(w_idx) != w_tail_outward)
							continue;
						const unsigned char mark = 
							vertex_marks[e2.get_other_v_idx(w_idx)];
						if ((mark & INPLAY)
						    && ((mark & INPLAY_OUTWARD)
							    ? e2.other_v_outward(w_idx) :
							      !e2.other_v_outward(w_idx)))
						{
							vertex_marks[e2.get_other_v_idx(w_idx)] =
									ELIMINATED;
						}
					}
				}
			}

			for (const edge_idx_t edge_idx : v.edge_indices())
				if (vertex_marks[edges[edge_idx].get_other_v_idx(v_idx)] == ELIMINATED)
					reduce_edge[edge_idx] = true;

			for (const edge_idx_t edge_idx : v.edge_indices())
				vertex_marks[edges[edge_idx].get_other_v_idx(v_idx)] = VACANT;
		cont:
			v_head_outward = !v_head_outward;
		} while (!v_head_outward);
	}

	info("Transitive reduction algorithm complete.  Now updating the string graph");

	// Update the directed string graph to remove the edges marked %true in
	// the @reduce_edge array.

	// Map from the old edge indices to the new edge indices.
	std::vector<edge_idx_t> new_edge_indices(edges.size());
	size_t i, j;
	for (i = 0, j = 0; i < edges.size(); i++) {
		if (reduce_edge[i]) {
			new_edge_indices[i] = std::numeric_limits<edge_idx_t>::max();
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

	// Re-number the edge indices list of each vertex, and remove any
	// indices that correspond to edges that were removed.
	for (BidirectedStringGraphVertex & v : vertices) {
		std::vector<edge_idx_t> & edge_indices = v.edge_indices();
		size_t j = 0;
		for (const edge_idx_t edge_idx : edge_indices) {
			const edge_idx_t new_edge_idx = new_edge_indices[edge_idx];
			if (new_edge_idx != std::numeric_limits<edge_idx_t>::max()) {
				edge_indices[j++] = new_edge_idx;
			}
		}
		edge_indices.resize(j);
	}

	info("Done removing transitive edges");
}
