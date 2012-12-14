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
// 	v >--> w  (v.E -> w.E, w.B -> v.E)
// 	v >--< w  (v.E -> w.B, w.E -> v.B)
// 	v <--< w  (v.B -> w.B, w.E -> v.E)
// 	v <--> w  (v.B -> w.E, w.B -> v.E)
void BidirectedStringGraph::transitive_reduction()
{
	unimplemented();

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
	for (v_idx_t v_idx = 0; v_idx < digraph.num_vertices(); v_idx++) {
		//
		// digraph edge
		// v.?1 -> w.?2   such that v_idx <= w_idx
		//      seq
		//
		// find edge w.^?1 -> v.^?2
		//
		// add bidigraph edge
		//
		//   v ?-? w where dirs =
		//   (?1 == E ? TAG_F_E : TAG_F_B) | (?2 == E ? TAG_G_E : TAG_G_B)
		//
		//   seq
		//
		const DirectedStringGraphVertex & v = digraph.vertices()[v_idx];
		foreach(const edge_idx_t v_w_edge_idx, v.edge_indices()) {
			const DirectedStringGraphEdge & v_w = digraph.edges()[v_w_edge_idx];
			const v_idx_t w_idx = v_w.get_v2_idx();

			/* Consider only 1 edge in each edge pair.
			 *
			 * - If an edge goes between different vertices, we can
			 *   skip edges where the first vertex has a higher
			 *   index than the second.
			 *
			 * - If an edge is a loop, then the edge pair is v.B ->
			 *   v.B and v.E -> v.E, so skip the v.B -> v.B edge. */
			if (v_idx < w_idx || (v_idx == w_idx && (v_idx & 1))) {
				const edge_idx_t w_v_edge_idx =
						digraph.locate_edge(w_idx ^ 1, v_idx ^ 1);

				const DirectedStringGraphEdge & w_v =
						digraph.edges()[w_v_edge_idx];

				v_idx_t dirs = ((v_idx & 1) << 1) | (w_idx & 1);

				BidirectedStringGraphEdge e;
				const v_idx_t v1_idx = v_idx / 2;
				const v_idx_t v2_idx = w_idx / 2;

				//e.get_seq_1_to_2().set_from_bv(v_g.get_seq());
				//e.get_seq_2_to_1().set_from_bv(w_f.get_seq());
				e.get_seq_1_to_2() = v_w.get_seq();
				e.get_seq_2_to_1() = w_v.get_seq();
				e.set_v_indices(v1_idx, v2_idx);
				e.set_dirs(dirs);
				e.set_mapped_read_count((w_v.get_mapped_read_count() +
							 v_w.get_mapped_read_count()) / 2);
				e.set_A_statistic((w_v.get_A_statistic() +
						   v_w.get_A_statistic()) / 2);
				e.set_traversal_count((w_v.get_traversal_count() +
						       v_w.get_traversal_count()));
				//if (w_v.get_num_inner_vertices() !=
				    //v_w.get_num_inner_vertices())
				//{
					//std::cerr << w_v << std::endl;
					//std::cerr << v_w << std::endl;
					//assert(0);
				//}
				e.set_num_inner_vertices((w_v.get_num_inner_vertices() +
							 v_w.get_num_inner_vertices()) / 2);
				if (v_w.is_special()) {
					assert(w_v.is_special());
					e.set_special();
				}

				edge_idx_t edge_idx = this->push_back_edge(e);
				_vertices[v1_idx].add_edge_idx(edge_idx);
				_vertices[v2_idx].add_edge_idx(edge_idx);
			}
		}
	}
	_orig_num_reads = digraph._orig_num_reads;
	info("Done building bidirected string graph from directed string graph");
}

void BidirectedStringGraph::print_stats(std::ostream & os) const
{
	os << "BidirectedStringGraph {" << std::endl;
	os << "    Number of vertices: " << num_vertices() << std::endl;
	os << "    Number of edges: " << num_edges() << std::endl;

	std::vector<unsigned char> out_degrees(num_vertices(), 0);
	std::vector<unsigned char> in_degrees(num_vertices(), 0);
	std::vector<v_idx_t> dir_histo(4, 0);
	foreach(const BidirectedStringGraphEdge & e, _edges) {
		v_idx_t v1_idx, v2_idx;
		e.get_v_indices(v1_idx, v2_idx);
		dir_histo[e.get_dirs()]++;
		if (e.v1_inward()) {
			if (in_degrees[v1_idx] + 1 != 0)
				in_degrees[v1_idx]++;
		} else {
			if (out_degrees[v1_idx] + 1 != 0)
				out_degrees[v1_idx]++;
		}
		if (e.v2_inward()) {
			if (in_degrees[v2_idx] + 1 != 0)
				in_degrees[v2_idx]++;
		} else {
			if (out_degrees[v2_idx] + 1 != 0)
				out_degrees[v2_idx]++;
		}
	}
	std::vector<v_idx_t> out_degree_hist(0xff, 0);
	std::vector<v_idx_t> in_degree_hist(0xff, 0);
	std::vector<v_idx_t> in_out_degree_hist(0xffff, 0);
	v_idx_t v_in_neq_out = 0;
	for (v_idx_t v_idx = 0; v_idx < num_vertices(); v_idx++) {
		out_degree_hist[out_degrees[v_idx]]++;
		in_degree_hist[in_degrees[v_idx]]++;
		if (out_degrees[v_idx] != in_degrees[v_idx])
			v_in_neq_out++;
		in_out_degree_hist[(v_idx_t(in_degrees[v_idx]) << 8) + out_degrees[v_idx]]++;
	}
	v_idx_t max_out_degree = 0, max_in_degree = 0;
	for (size_t i = 0; i < out_degree_hist.size(); i++)
		if (out_degree_hist[i] != 0)
			max_out_degree = i;
	for (size_t i = 0; i < in_degree_hist.size(); i++)
		if (in_degree_hist[i] != 0)
			max_in_degree = i;
	os << "    Number of isolated vertices: "
	   << in_out_degree_hist[0x000] << std::endl;
	os << "    Number of inner vertices: "
	   << in_out_degree_hist[0x101] << std::endl;
	os << "    Number of branch beginning vertices: "
	   << in_out_degree_hist[0x001] << std::endl;
	os << "    Number of branch ending vertices: "
	   << in_out_degree_hist[0x100] << std::endl;
	os << "    Number of vertices with unequal in degree and out degree: "
	   << v_in_neq_out << std::endl;
	os << "    Max in degree: " << max_in_degree
	   << (max_in_degree == 0xff ? '+' : ' ') << std::endl;
	os << "    Max out degree: " << max_out_degree
	   << (max_out_degree == 0xff ? '+' : ' ') << std::endl;

	// 11
	// (v1_outward, v2_inward)
	// 00 => in, out    <---<
	// 01 => in, in     <--->
	// 10 => out, out   >---<
	// 11 => out, in    >--->
	os << "    Number of edges >--->: " << (dir_histo[0x0] + dir_histo[0x3]) << std::endl;
	os << "    Number of edges <--->: " << (dir_histo[0x1]) << std::endl;
	os << "    Number of edges >---<: " << (dir_histo[0x2]) << std::endl;
	os << "}" << std::endl;
}

void BidirectedStringGraph::map_contained_read(const v_idx_t downstream_read_idx,
					       const v_idx_t downstream_read_dir,
					       const BaseVec::size_type overhang_len)
{
	unimplemented();
}

void BidirectedStringGraph::calculate_A_statistics()
{
	unimplemented();
}

void BidirectedStringGraph::min_cost_circulation()
{
	unimplemented();
}

void BidirectedStringGraph::assert_eulerian_cycle_valid(const std::vector<size_t> & cycle) const
{
	{
		std::vector<edge_idx_t> times_traversed(num_edges(), 0);
		foreach (const size_t edge_idx, cycle) {
			times_traversed[edge_idx]++;
		}
		for (size_t i = 0; i < num_edges(); i++) {
			assert(times_traversed[i] == _edges[i].get_traversal_count());
		}
	}
	if (cycle.empty())
		return;

	const BidirectedStringGraphEdge & first_edge = _edges[cycle[0]];
	v_idx_t cur_v_idx = first_edge.get_v1_idx();

	// Direction of head through which cur_v_idx was entered
	bool inward = !first_edge.v1_inward();

	if (cycle.size() > 1) {
		v_idx_t second_v_idx = first_edge.get_v2_idx();
		// We want to go around the cycle in the direction indicated.
		// So, for the first edge between vertices A and B, check to
		// make sure the second edge in the cycle includes vertex B with
		// the direction of the adjacent head consistent with the first
		// head.  If not, flip vertices B and A so we go the other way
		// around the cycle.
		const BidirectedStringGraphEdge & second_edge = _edges[cycle[1]];
		if ((second_edge.get_v1_idx() != second_v_idx ||
		     second_edge.v1_inward() == first_edge.v2_inward())
		    && (second_edge.get_v2_idx() != second_v_idx ||
			second_edge.v2_inward() == first_edge.v2_inward()))
		{
			cur_v_idx = first_edge.get_v2_idx();
			inward = !first_edge.v2_inward();
		}
	}

	for (size_t i = 0; i < cycle.size(); i++) {
		// Given that cur_v_idx has been entered through a head in
		// direction @inward relative to @cur_v_idx, assert that the
		// next edge e in the cycle contains vertex cur_v_idx and can be
		// traversed.
		info("cur_v_idx = %u", cur_v_idx);
		const BidirectedStringGraphEdge & e = _edges[i];

		e.print(std::cout, cur_v_idx);

		assert(cur_v_idx == e.get_v1_idx() ||
		       cur_v_idx == e.get_v2_idx());

		if (e.get_v1_idx() == cur_v_idx &&
		    e.v1_inward() != inward)
		{
			inward = e.v2_inward();
			cur_v_idx = e.get_v2_idx();
			continue;
		} 
		if (e.get_v2_idx() == cur_v_idx &&
		    e.v2_inward() != inward)
		{
			inward = e.v1_inward();
			cur_v_idx = e.get_v1_idx();
			continue;
		}
		assert(0);
	}
	assert(inward != _edges[cycle[0]].v2_inward());
	assert(cur_v_idx == _edges[cycle[0]].get_v2_idx());
}

// XXX this only checks if in-degree == out-degree for all vertices.  The graph
// must also be strongly connected for an eulerian cycle to exist.
void BidirectedStringGraph::assert_eulerian_cycle_possible() const
{
	for (v_idx_t v_idx = 0; v_idx < num_vertices(); v_idx++) {
		const BidirectedStringGraphVertex & v = _vertices[v_idx];
		unsigned long in_degree = 0;
		unsigned long out_degree = 0;
		foreach (const edge_idx_t edge_idx, v.edge_indices()) {
			const BidirectedStringGraphEdge & e = _edges[edge_idx];
			v_idx_t v1_idx, v2_idx;
			e.get_v_indices(v1_idx, v2_idx);
			if (v_idx == v1_idx) {
				if (e.v1_inward())
					in_degree += e.get_traversal_count();
				else
					out_degree += e.get_traversal_count();
			}
			if (v_idx == v2_idx) {
				if (e.v2_inward())
					in_degree += e.get_traversal_count();
				else
					out_degree += e.get_traversal_count();
			}
		}
		assert(in_degree == out_degree);
	}
}

void BidirectedStringGraph::eulerian_cycle(std::vector<size_t> & cycle) const
{
	assert_eulerian_cycle_possible();

	v_idx_t n_verts = num_vertices();
	edge_idx_t n_edges = num_edges();

	info("Finding Eulerian path in bidirected graph");
	info("num_vertices = %lu", n_verts);
	info("num_edges = %lu", n_edges);

	unsigned long total_traversal_count = 0;
	unsigned long num_special_edges = 0;
	unsigned long special_traversal_count = 0;
	v_idx_t start_v_idx;
	bool found_start_v_idx = false;
	foreach (const BidirectedStringGraphEdge & e, _edges) {
		if (e.is_special()) {
			num_special_edges++;
			special_traversal_count += e.get_traversal_count();
		}
		total_traversal_count += e.get_traversal_count();
		if (!found_start_v_idx) {
			if (e.v1_outward()) {
				start_v_idx = e.get_v1_idx();
				found_start_v_idx = true;
			} else if (e.v2_outward()) {
				start_v_idx = e.get_v2_idx();
				found_start_v_idx = true;
			}
		}
	}
	info("total_traversal_count = %lu", total_traversal_count);
	info("num_special_edges = %lu", num_special_edges);
	info("special_traversal_count = %lu", special_traversal_count);

	if (!found_start_v_idx) {
		info("WARNING: Empty Eulerian path!");
		return;
	}
	info("Starting at vertex %lu", start_v_idx + 1);

	// Start with empty path and empty stack
	cycle.clear();
	cycle.reserve(total_traversal_count);

	struct stack_elem {
		// index of vertex
		v_idx_t v_idx;

		// index of edge we're taking out of vertex v_idx
		edge_idx_t edge_idx;

		// direction of head relative to vertex v_idx through which
		// vertex v_idx was entered
		bool inward;
	};
	std::vector<stack_elem> stack;

	std::vector<edge_idx_t> in_indices(n_verts, 0);
	std::vector<edge_idx_t> out_indices(n_verts, 0);
	std::vector<unsigned> times_traversed(n_edges, 0);

	// Start by looking for out heads; this implies that we consider the
	// starting vertex to have been entered through an inward head.
	bool inward = true;
	v_idx_t v_idx = start_v_idx;
	while (1) {
		edge_idx_t edge_idx;
		bool pop_stack = true;
		const BidirectedStringGraphVertex & v = _vertices[v_idx];
		const BidirectedStringGraphEdge *e;
		if (inward) {
			// v was entered through an inwards head, so look for
			// adjacent edges that still have traversal count
			// remaining and begin with a head pointed outwards
			// relative to v.
			for (; out_indices[v_idx] < v.degree(); out_indices[v_idx]++) {
				edge_idx = v.edge_indices()[out_indices[v_idx]];
				e = &_edges[edge_idx];

				if (times_traversed[edge_idx] < e->get_traversal_count() &&
				    e->v_outward(v_idx))
				{
					pop_stack = false;
					break;
				}
			}
		} else {
			// v was entered through an outwards head, so look for
			// adjacent edges that still have traversal count
			// remaining and begin with a head pointed inwards
			// relative to v.
			for (; in_indices[v_idx] < v.degree(); in_indices[v_idx]++) {
				edge_idx = v.edge_indices()[in_indices[v_idx]];
				e = &_edges[edge_idx];
				if (times_traversed[edge_idx] < e->get_traversal_count() &&
				    e->v_inward(v_idx))
				{
					pop_stack = false;
					break;
				}
			}
		}
		if (pop_stack) {

			// v does not have an edge that can be traversed.
			//
			// Pop (v_idx, edge_idx, inward) from stack and add
			// edge_idx to cycle.  If stack is empty, the cycle is
			// complete.
			if (stack.empty())
				break;
			const stack_elem & elem = stack.back();
			v_idx = elem.v_idx;
			edge_idx = elem.edge_idx;
			inward = elem.inward;
			stack.pop_back();

			info("Popped (v_idx=%u, edge_idx=%u, inward=%d)",
			     v_idx, edge_idx, inward);
			_edges[edge_idx].print(std::cout, v_idx, false);

			cycle.push_back(edge_idx);
		} else {
			// v has an edge (_edges[edge_idx]) that can be traversed.
			//
			// Push (v_idx, edge_idx, inward) onto the stack, then
			// mark the edge as traversed one additional time, then
			// move to the vertex at the other end of the edge by
			// updating v_idx and edge_idx.
			stack_elem elem;
			elem.v_idx = v_idx;
			elem.edge_idx = edge_idx;
			elem.inward = inward;

			info("Pushed (v_idx=%u, edge_idx=%u, inward=%d)",
			     v_idx, edge_idx, inward);
			stack.push_back(elem);
			times_traversed[edge_idx]++;

			v_idx = _edges[edge_idx].get_other_v_idx(v_idx);
			if (_edges[edge_idx].is_loop()) {
				inward ^= (_edges[edge_idx].v1_inward() ==
					   _edges[edge_idx].v2_inward());
			} else {
				inward = _edges[edge_idx].this_v_inward(v_idx);
			}
		}
	}
	assert_eulerian_cycle_valid(cycle);
}
