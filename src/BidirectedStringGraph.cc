#include "BidirectedStringGraph.h"
#include "DirectedStringGraph.h"

const char BidirectedStringGraph::magic[]
		= {'B', 'i', 'd', 'i', 'g', 'r', 'a', 'p', 'h', '\0'};

//
// Transitive reduction of a bidirected string graph.
//
// This is based on the code for transitive reduction of a directed string graph
// (DirectedStringGraph::transitive_reduction()), but some modifications need to
// be made to make it work on bidirected graphs.
//
// In particular, an edge v ?-? x can only be considered transitive given a pair
// of edges v ?-? w ?-? x if:
// 	(1) The two heads adjacent to w have opposite orientation;
// 	(2) The heads adjacent to v in v ?-? w and v ?-? x have the same
// 	orientation, and
// 	(3) The heads adjacent to x in v ?-? x and w ?-? x have the same
// 	orientation.
// 	(4) And as in the directed case, the sequence v ?-? x must be equal to
// 	the sequence v ?-? w + v ?-? x.
//
// In addition, in the directed case, there can be at most one edge v -> w
// (although there may be an edge w -> v as well).  But in the bidirected case,
// there can be up to 4 edges between v and w:
// 	v >--> w
// 	v <--< w
// 	v >--< w
// 	v <--> w
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
		// First, only consider edges leaving v with an out head.  Then,
		// on the second iteration, only consider edges leaving v with
		// an in head.
		do {
			const BidirectedStringGraphVertex & v = vertices[v_idx];

			BaseVec::size_type longest = 0;

			foreach(const edge_idx_t edge_idx, v.edge_indices()) {
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

			foreach(const edge_idx_t edge_idx, v.edge_indices()) {
				const v_idx_t w_idx = edges[edge_idx].get_other_v_idx(v_idx);
				if (vertex_marks[w_idx] & INPLAY) {
					const bool w_tail_outward =
							edges[edge_idx].other_v_outward(v_idx);
					const BidirectedStringGraphVertex & w = vertices[w_idx];
					foreach(const edge_idx_t w_edge_idx, w.edge_indices()) {
						const BidirectedStringGraphEdge & e2 = edges[w_edge_idx];
						if (e2.length() > longest)
							break;
						if (e2.this_v_outward(w_idx) == w_tail_outward)
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

			foreach(const edge_idx_t edge_idx, v.edge_indices()) {
				const BidirectedStringGraphEdge & e = edges[edge_idx];
				if (e.this_v_outward(v_idx) == v_head_outward
				    && vertex_marks[e.get_other_v_idx(v_idx)] == ELIMINATED)
				{
						reduce_edge[edge_idx] = true;
				}
			}

			foreach(const edge_idx_t edge_idx, v.edge_indices())
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
	foreach(BidirectedStringGraphVertex & v, vertices) {
		std::vector<edge_idx_t> & edge_indices = v.edge_indices();
		size_t j = 0;
		foreach(const edge_idx_t edge_idx, edge_indices) {
			const edge_idx_t new_edge_idx = new_edge_indices[edge_idx];
			if (new_edge_idx != std::numeric_limits<edge_idx_t>::max()) {
				edge_indices[j++] = new_edge_idx;
			}
		}
		edge_indices.resize(j);
	}

	info("Done removing transitive edges");
}

void BidirectedStringGraph::collapse_unbranched_paths()
{
	unimplemented();
}

void BidirectedStringGraph::build_from_digraph(const DirectedStringGraph & digraph)
{
	info("Building bidirected string graph from directed string graph");
	info("Directed string graph: %zu vertices", digraph.num_vertices());
	info("Bidirected string graph: %zu vertices", this->num_vertices());
	assert(digraph.num_vertices() % 2 == 0);
	assert(num_vertices() == digraph.num_vertices() / 2);
	assert(num_edges() == 0);
	for (v_idx_t f_idx = 0; f_idx < digraph.num_vertices(); f_idx++) {
		//
		// digraph edge
		// f.?1 -> g.?2   such that f_idx <= g_idx
		//      seq
		//
		// find edge g.^?1 -> f.^?2
		// 
		// add bidigraph edge
		//
		//   f ?-? g where dirs =
		//   (?1 == E ? TAG_F_E : TAG_F_B) | (?2 == E ? TAG_G_E : TAG_G_B)
		//
		//   seq
		//
		foreach(const edge_idx_t f_g_edge_idx, digraph.vertices()[f_idx].edge_indices())
		{
			const DirectedStringGraphEdge & f_g = digraph.edges()[f_g_edge_idx];
			const v_idx_t g_idx = f_g.get_v2_idx();

			if (f_idx == g_idx)
				unimplemented();
			if (f_idx < g_idx) {
				const edge_idx_t g_f_edge_idx =
						digraph.locate_edge(g_idx ^ 1, f_idx ^ 1);

				const DirectedStringGraphEdge & g_f =
						digraph.edges()[g_f_edge_idx];

				v_idx_t dirs = 0;
				dirs |= (f_idx & 1) ?  TAG_F_E : TAG_F_B;
				dirs |= (g_idx & 1) ?  TAG_G_E : TAG_G_B;

				BidirectedStringGraphEdge e;
				const v_idx_t v1_idx = f_idx / 2;
				const v_idx_t v2_idx = g_idx / 2;

				//e.get_seq_1_to_2().set_from_bv(f_g.get_seq());
				//e.get_seq_2_to_1().set_from_bv(g_f.get_seq());
				e.get_seq_1_to_2() = f_g.get_seq();
				e.get_seq_2_to_1() = g_f.get_seq();
				e.set_v_indices(v1_idx, v2_idx);
				e.set_dirs(dirs);

				edge_idx_t edge_idx = this->push_back_edge(e);
				_vertices[v1_idx].add_edge_idx(edge_idx);
				_vertices[v2_idx].add_edge_idx(edge_idx);
			}
		}
	}
	info("Done building bidirected string graph from directed string graph");
}
