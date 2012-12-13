#include "DirectedStringGraph.h"
#include <fstream>
#include <ostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "BidirectedStringGraph.h"
#include <math.h>
#include <lemon/network_simplex.h>
#include <lemon/smart_graph.h>

const char DirectedStringGraph::magic[] =
	{'D', 'i', 'g', 'r', 'a', 'p', 'h', '\0', '\0', '\0'};

//
// Perform transitive edge reduction on a directed string graph:
// Remove all edges v -> x where there exist edges v -> w -> x.
//
// See:
//    "The fragment string assembly graph", Eugene W. Myers, 2005
//
// There are, however, some changes that have to be made to make the algorithm
// work in practice.  In particular, the edge v -> x should only be removed if
// it is actually labeled by the same sequence as v -> w -> x.
//
// It's important to keep in mind that there can only at most one edge v -> w
// for any v and w, since this represents one kind of overlap.  But, it's
// possible for there to be an edge w -> v as well.
//
void DirectedStringGraph::transitive_reduction()
{
	info("Performing transitive reduction on directed string graph with "
	     "%zu vertices and %zu edges", this->num_vertices(), this->num_edges());

	std::vector<DirectedStringGraphVertex> & vertices = this->vertices();
	std::vector<DirectedStringGraphEdge> & edges = this->edges();

	info("Sorting adjacency lists of vertices by edge length");
	this->sort_adjlists_by_edge_len();

	static const unsigned char VACANT = 0;
	static const unsigned char INPLAY = 1;
	static const unsigned char ELIMINATED = 2;

	// Per-vertex mark.  Initially set to VACANT.
	std::vector<unsigned char> vertex_marks(this->num_vertices(), VACANT);

	// Map from vertex index to index of edges back to the vertex currently
	// under consideration
	std::vector<edge_idx_t> v_idx_to_back_edge_idx(this->num_vertices(),
						       std::numeric_limits<edge_idx_t>::max());

	// Per-edge boolean indicating whether each edge is to be deleted as
	// part of the transitive reduction or not.
	std::vector<bool> reduce_edge(this->num_edges(), false);

	info("Looking for transitive edges based at each of %zu vertices",
	     vertices.size());

	// Iterate through every vertex @v in the graph that has outgoing edges.
	for (size_t v_idx = 0; v_idx < vertices.size(); v_idx++) {
		const DirectedStringGraphVertex & v = vertices[v_idx];

		if (v.out_degree() == 0)
			continue;

		// Mark each vertex adjacent to @v as INPLAY, and initialize the
		// map from the adjacent vertices' indices to the back edges'
		// indices.
		foreach(const edge_idx_t edge_idx, v.edge_indices()) {
			const DirectedStringGraphEdge & e = edges[edge_idx];
			const v_idx_t w_idx = e.get_v2_idx();
			assert(edge_idx != std::numeric_limits<edge_idx_t>::max());
			assert(v_idx_to_back_edge_idx[w_idx] ==
			       std::numeric_limits<edge_idx_t>::max());
			v_idx_to_back_edge_idx[w_idx] = edge_idx;

			vertex_marks[w_idx] = INPLAY;
		}

		// Length of the longest sequence label on the edges leaving
		// vertex @v.
		const BaseVec::size_type longest = edges[v.edge_indices().back()].length();

		// For each outgoing edge from v -> w in order of labeled
		// sequence length, consider each vertex w that is still marked
		// INPLAY.
		foreach(const edge_idx_t edge_idx, v.edge_indices()) {
			const DirectedStringGraphEdge & e = edges[edge_idx];
			const v_idx_t w_idx = e.get_v2_idx();

			if (vertex_marks[w_idx] != INPLAY)
				continue;

			// The edge v -> w must be an irreducible edge if w is
			// still marked INPLAY at this point, since all shorter
			// edges were already considered.
			//
			// Now, consider the edges leaving vertex w.  Each such
			// edge that goes to a vertex marked INPLAY must be
			// directly reachable from v, and therefore the edge
			// must be removed, unless its sequence does not
			// actually match the sequence from the edges
			// v -> w -> x.
			const DirectedStringGraphVertex & w = vertices[w_idx];
			foreach(const edge_idx_t w_edge_idx, w.edge_indices()) {
				const DirectedStringGraphEdge & e2 = edges[w_edge_idx];
				if (e.length() + e2.length() > longest)
					break;

				const v_idx_t x_idx = e2.get_v2_idx();
				if (vertex_marks[x_idx] != INPLAY)
					continue;

				const edge_idx_t back_edge_idx =
						v_idx_to_back_edge_idx[x_idx];

				assert(back_edge_idx !=
				       std::numeric_limits<edge_idx_t>::max());

				const DirectedStringGraphEdge & back_edge =
					edges[back_edge_idx];

				assert(back_edge.get_v1_idx() == v_idx);
				assert(back_edge.get_v2_idx() == x_idx);


				if (e.length() + e2.length() != back_edge.length())
					continue;

				for (BaseVec::size_type i = 0; i < e.length(); i++)
					if (e.get_seq()[i] != back_edge.get_seq()[i])
						goto next_edge;

				for (BaseVec::size_type i = 0; i < e2.length(); i++)
					if (e2.get_seq()[i] !=
					    back_edge.get_seq()[i + e.length()])
						goto next_edge;

				vertex_marks[x_idx] = ELIMINATED;
				next_edge:
				;
			}
		}

		// Once again, go through the outgoing edges from v.  For each
		// neighboring vertex marked ELIMINATED, mark the corresponding
		// edge(s) for reduction.  Return both INPLAY and ELIMINATED
		// vertices to VACANT status.
		foreach(const edge_idx_t edge_idx, v.edge_indices()) {
			const DirectedStringGraphEdge & e = edges[edge_idx];
			const v_idx_t w_idx = e.get_v2_idx();
			if (vertex_marks[w_idx] == ELIMINATED)
				reduce_edge[edge_idx] = true;
			v_idx_to_back_edge_idx[w_idx] = std::numeric_limits<edge_idx_t>::max();
			vertex_marks[w_idx] = VACANT;
		}
	}

	info("Transitive reduction algorithm complete.  Now updating the string graph");

	for (size_t i = 0; i < edges.size(); i++) {
		if (reduce_edge[i]) {
			const DirectedStringGraphEdge & e = _edges[i];
			const edge_idx_t j = locate_edge(e.get_v2_idx() ^ 1,
							 e.get_v1_idx() ^ 1);
			if (!reduce_edge[j]) {
				std::cout << "The following 2 edges are opposites "
					<< "but were not both reduced:" << std::endl;
				e.print(std::cout, 0, true);
				std::cout << std::endl;
				_edges[j].print(std::cout, 0, true);
				std::cout << std::endl;
			}
		}
	}

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
	const size_t num_remaining_edges = j;
	edges.resize(num_remaining_edges);

	info("Removing %zu of %zu edges (%.2f%%)",
	     num_removed_edges, num_original_edges,
	     TO_PERCENT(num_removed_edges, num_original_edges));

	// Re-number the edge indices list of each vertex, and remove any
	// indices that correspond to edges that were removed.
	foreach(DirectedStringGraphVertex & v, vertices) {
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

void DirectedStringGraph::follow_unbranched_path(DirectedStringGraphEdge & e,
						 std::vector<bool> & remove_edge,
						 std::vector<bool> & remove_vertex,
						 const std::vector<bool> & v_inner)
{
	BaseVecVec::size_type new_seq_len = e.length();
	v_idx_t vi_idx = e.get_v2_idx();
	assert(v_inner[vi_idx]);

	v_idx_t num_inner_vertices = 0;
	// Found beginning of unbranched path.  Walk along it until
	// the end to get the total sequence length.
	do {
		const DirectedStringGraphVertex &vi = _vertices[vi_idx];
		assert(vi.out_degree() == 1);
		const DirectedStringGraphEdge &ei_i1 = _edges[vi.first_edge_idx()];
		if (new_seq_len + ei_i1.length() < new_seq_len)
			fatal_error("Edge too long");
		new_seq_len += ei_i1.length();
		num_inner_vertices++;
		vi_idx = ei_i1.get_v2_idx();
	} while (v_inner[vi_idx]);

	// Now build the new sequence, update the first edge in the path to skip
	// to the last vertex in the path, and mark the other edges for removal
	// in the @remove_edge array.
	BaseVec & new_seq = e.get_seq();
	BaseVec::size_type seq_idx = new_seq.length();
	new_seq.resize(new_seq_len);
	vi_idx = e.get_v2_idx();

	e.set_num_inner_vertices(num_inner_vertices);
	while (num_inner_vertices--) {
		const DirectedStringGraphVertex &vi = _vertices[vi_idx];
		const edge_idx_t ei_i1_idx = vi.first_edge_idx();
		const DirectedStringGraphEdge &ei_i1 = _edges[ei_i1_idx];
		const BaseVec & ei_i1_seq = ei_i1.get_seq();
		for (BaseVec::size_type i = 0; i < ei_i1_seq.length(); i++) {
			assert2(seq_idx < new_seq_len);
			new_seq.set(seq_idx++, ei_i1_seq[i]);
		}
		e.increment_mapped_read_count(ei_i1.get_mapped_read_count());
		remove_edge[ei_i1_idx] = true;
		remove_vertex[vi_idx] = true;
		vi_idx = ei_i1.get_v2_idx();
	}
	assert(!v_inner[vi_idx]);
	assert(new_seq_len == seq_idx);
	e.set_v2_idx(vi_idx);
}

void DirectedStringGraph::collapse_unbranched_paths()
{
	const v_idx_t n_verts = num_vertices();
	const edge_idx_t n_edges = num_edges();

	info("Collapsing unbranched paths in directed string graph");
	info("Original graph has %zu vertices and %zu edges", n_verts, n_edges);

	size_t num_inner_vertices = 0;
	// Find whether each vertex is inner or not.  A vertex is inner iff it
	// has indegree 1 and outdegree 1.
	std::vector<bool> v_inner(n_verts, false);
	{
		std::vector<unsigned char> v_in_degrees(n_verts, 0);
		std::vector<unsigned char> v_out_degrees(n_verts, 0);

		foreach (const DirectedStringGraphEdge & e, _edges) {
			v_idx_t v1_idx, v2_idx;
			e.get_v_indices(v1_idx, v2_idx);
			if (v_out_degrees[v1_idx] < 2)
				v_out_degrees[v1_idx]++;
			if (v_in_degrees[v2_idx] < 2)
				v_in_degrees[v2_idx]++;
		}
		for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
			if (v_in_degrees[v_idx] == 1 && v_out_degrees[v_idx] == 1) {
				v_inner[v_idx] = true;
				num_inner_vertices++;
			}
		}
	}

	info("Found %zu inner vertices (%.2f%% of all vertices)",
	     num_inner_vertices, TO_PERCENT(num_inner_vertices, n_verts));

	// Go through each non-inner vertex and look for any neighboring inner
	// vertices.  These are the starts of unbranched paths that will be
	// collapsed.  For each such path, follow it until its end and determine
	// the total length of the edges.  Then, set the first edge's label to
	// the concatenation of the labels of all the edges in the unbranched
	// path, then make the first edge point directly to the last vertex in
	// the path and mark the other edges for deletion.
	//
	// Note: smooth rings are not collapsed yet.
	//
	size_t num_unbranched_paths = 0;
	std::vector<bool> remove_edge(n_edges, false);
	std::vector<bool> remove_vertex(n_verts, false);
	for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
		if (!v_inner[v_idx]) {
			const DirectedStringGraphVertex & v = _vertices[v_idx];
			foreach(edge_idx_t edge_idx, v.edge_indices()) {
				DirectedStringGraphEdge & e = _edges[edge_idx];
				v_idx_t v2_idx = e.get_v2_idx();
				if (v_inner[v2_idx] && !remove_vertex[v2_idx]) {
					num_unbranched_paths++;
					follow_unbranched_path(e, remove_edge,
							       remove_vertex, v_inner);
				}
			}
		}
	}

	for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
		if (v_inner[v_idx] && !remove_vertex[v_idx]) {
			std::cerr << "Graph contains a smooth ring!" << std::endl;
			unimplemented();
		}
	}

	info("Found %zu unbranched paths", num_unbranched_paths);

	// Compute the new vertex indices.
	std::vector<v_idx_t> old_to_new_v_indices(n_verts,
						  std::numeric_limits<v_idx_t>::max());
	v_idx_t new_v_idx = 0;
	for (v_idx_t old_v_idx = 0; old_v_idx < n_verts; old_v_idx++)
		if (!remove_vertex[old_v_idx])
			old_to_new_v_indices[old_v_idx] = new_v_idx++;

	info("Updated vertices are indexed [0, %lu)", new_v_idx);

	info("Updating edges");
	// Compute the new edge indices, set the new vertex indices in each edge,
	// and move the edges
	std::vector<edge_idx_t> old_to_new_edge_indices(n_edges,
							std::numeric_limits<edge_idx_t>::max());
	edge_idx_t new_edge_idx = 0;
	for (edge_idx_t old_edge_idx = 0; old_edge_idx < n_edges; old_edge_idx++) {
		if (!remove_edge[old_edge_idx]) {
			DirectedStringGraphEdge & e = _edges[old_edge_idx];
			old_to_new_edge_indices[old_edge_idx] = new_edge_idx;
			v_idx_t v1_idx, v2_idx;
			e.get_v_indices(v1_idx, v2_idx);
			e.set_v_indices(old_to_new_v_indices[v1_idx],
					old_to_new_v_indices[v2_idx]);
			assert2(e.get_v1_idx() != std::numeric_limits<v_idx_t>::max() &&
				e.get_v2_idx() != std::numeric_limits<v_idx_t>::max());
			_edges[new_edge_idx++] = _edges[old_edge_idx];
		}
	}
	info("Updated edges are indexed [0, %lu)", new_edge_idx);
	info("%zu edges were removed (%f%% of total)",
	     _edges.size() - new_edge_idx,
	     TO_PERCENT(_edges.size() - new_edge_idx, _edges.size()));

	_edges.resize(new_edge_idx);

	info("Updating vertices");
	// Set new edge indices in each vertex and move the vertices
	new_v_idx = 0;
	for (v_idx_t old_v_idx = 0; old_v_idx < n_verts; old_v_idx++) {
		if (!v_inner[old_v_idx]) {
			DirectedStringGraphVertex & v = _vertices[old_v_idx];
			std::vector<edge_idx_t> & edge_indices = v.edge_indices();
			for (size_t i = 0; i < v.out_degree(); i++) {
				edge_indices[i] = old_to_new_edge_indices[edge_indices[i]];
				assert2(edge_indices[i] != std::numeric_limits<edge_idx_t>::max());
			}
			_vertices[new_v_idx++] = _vertices[old_v_idx];
		}
	}
	assert(new_v_idx == n_verts - num_inner_vertices);
	_vertices.resize(new_v_idx);
	info("Done collapsing unbranched paths in directed string graph");
}

void DirectedStringGraph::mark_component(const v_idx_t v_idx,
					 std::vector<bool> & visited,
					 v_idx_t & component_size) const
{
	assert(!visited[v_idx]);
	visited[v_idx] = true;
	component_size++;
	foreach (edge_idx_t edge_idx, _vertices[v_idx].edge_indices()) {
		const v_idx_t w_idx = _edges[edge_idx].get_v2_idx();
		if (!visited[w_idx]) {
			mark_component(w_idx, visited, component_size);
		}
	}
}

void DirectedStringGraph::print_stats(std::ostream & os) const
{
	os << "DirectedStringGraph {" << std::endl;
	os << "    Number of vertices: " << num_vertices() << std::endl;
	os << "    Number of edges: " << num_edges() << std::endl;

	std::vector<unsigned char> out_degrees(num_vertices(), 0);
	std::vector<unsigned char> in_degrees(num_vertices(), 0);
	double total_mapped_count = 0.0;
	size_t more_than_one_mapped_count = 0;
	foreach(const DirectedStringGraphEdge & e, _edges) {
		if (out_degrees[e.get_v1_idx()] + 1 != 0) {
			out_degrees[e.get_v1_idx()]++;
		}
		if (in_degrees[e.get_v2_idx()] + 1 != 0) {
			in_degrees[e.get_v2_idx()]++;
		}
		assert(e.get_mapped_read_count() >= 1.0);
		total_mapped_count += e.get_mapped_read_count();
		if (e.get_mapped_read_count() > 1.0)
			more_than_one_mapped_count++;
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
	   << (max_in_degree == 0xff ? '+' : ' ') << std::endl;

	os << "    Number of edges that one or more contained reads map onto: "
	   << more_than_one_mapped_count << std::endl;
	os << "    Average number of contained reads that map onto each edge: "
	   << (num_edges() ? total_mapped_count / num_edges() : 0) << std::endl;

	std::vector<bool> visited(num_vertices(), false);
	std::vector<v_idx_t> component_sizes;

	for (v_idx_t v_idx = 0; v_idx < num_vertices(); v_idx++) {
		if (!visited[v_idx] && in_degrees[v_idx] == 0) {
			v_idx_t component_size = 0;
			mark_component(v_idx, visited, component_size);
			component_sizes.push_back(component_size);
		}
	}
	for (v_idx_t v_idx = 0; v_idx < num_vertices(); v_idx++) {
		if (!visited[v_idx]) {
			v_idx_t component_size = 0;
			mark_component(v_idx, visited, component_size);
			component_sizes.push_back(component_size);
		}
	}
	std::sort(component_sizes.begin(), component_sizes.end());
	os << "    Number of components: " << component_sizes.size() << std::endl;
	os << "    Component sizes:" << std::endl;
	foreach (v_idx_t component_size, component_sizes) {
		os << "        " << component_size << ' '
		   << (component_size == 1 ? "vertex" : "vertices") << std::endl;
	}
	os << "}" << std::endl;
}

//
// Consider a read f contained in another read g:
//
//
//                  f.B ---------> f.E
//          g.B ------------------------> g.E
//
// where the relevant portion of the actual graph might look like:
//
// 	v1.? -> v2.? -> v3.? -> g.E
//
// Or:
//
//                  f.B ---------> f.E
//          g.E <------------------------ g.B
//
//  where the relevant portion of the actual graph might look like:
//
//            v1.? -> v2.? -> v3.? -> g.B
//
// In reality, there are no vertices f.B and f.E because they were removed from
// the graph.  The goal here is to start at vertex g.E (or g.B) and walk
// backwards (e.g. towards the v?.? shown above) to find all possible edges into
// which the f.E vertex would be located, so that the number of reads that map
// onto those edges can be incremented.  Ideally, there would be only 1 such
// edge--- but due to branching, there may be many.  At most 100 mapped edges
// are considered before throwing away this contained read.
//
// @v_idx:
// 	Next branchpoint vertex.  (initially: g.E for forward case; g.B for
// 	reverse-complement case).
//
// @overhang_len:
// 	Number of bp of uncontained read remaining after threading the contained
// 	read all the way to the vertex of index @v_idx.
//
// @mapped_edges:
// 	Array of indices of edges onto which the contained read would map.
//
// @num_mapped_edges:
// 	Number of filled members of @mapped_edges.
//
// @max_mapped_edges:
// 	Size of max_mapped_edges.  The search is terminated when
// 	@num_mapped_edges reaches @max_mapped_edges.
//
void DirectedStringGraph::walk_back_edges(const v_idx_t v_idx,
					  size_t overhang_len,
					  edge_idx_t mapped_edges[],
					  size_t & num_mapped_edges,
					  size_t max_mapped_edges)
{
	assert(_back_edges.size() == _vertices.size());

	// Map the read onto any edges that go into the current vertex with
	// length less than the remaining read length.
	foreach (edge_idx_t edge_idx, _back_edges[v_idx]) {
		const DirectedStringGraphEdge & e = _edges[edge_idx];
		assert(e.get_v2_idx() == v_idx);
		if (overhang_len < e.length()) {
			mapped_edges[num_mapped_edges++] = edge_idx;
			if (num_mapped_edges == max_mapped_edges)
				return;
		}
	}

	// For each edge going into the current vertex with length greater than
	// or equal to the remaining read length, recursively call this
	// procedure on them with the remaining read length decremented by the
	// legnth of the edge.
	//
	// If upon return of the recursive call, the number of mapped edges has
	// reached the maximum, return early.
	foreach (edge_idx_t edge_idx, _back_edges[v_idx]) {
		const DirectedStringGraphEdge & e = _edges[edge_idx];
		if (overhang_len >= e.length()) {
			walk_back_edges(e.get_v1_idx(),
					overhang_len - e.length(),
					mapped_edges,
					num_mapped_edges,
					max_mapped_edges);
			if (num_mapped_edges == max_mapped_edges)
				return;
		}
	}
}

//
//
void DirectedStringGraph::map_contained_read(const v_idx_t downstream_read_idx,
					     const v_idx_t downstream_read_dir,
					     const BaseVec::size_type overhang_len)
{
	assert(downstream_read_idx < num_vertices() / 2);
	assert(downstream_read_dir < 2);

	const v_idx_t downstream_v_idx = downstream_read_idx * 2 + downstream_read_dir;

	// We need to be able to walk the graph in the opposite direction that
	// the edges are going--- so index the back edges, so that it's possible
	// to get a list of edges entering a vertex, given the vertex index.
	if (_back_edges.size() == 0) {
		info("Indexing back edges (num_vertices = %zu, num_edges = %zu)",
				num_vertices(), num_edges());
		_back_edges.resize(num_vertices());
		foreach (const DirectedStringGraphVertex & v, _vertices)
			foreach (const edge_idx_t edge_idx, v.edge_indices())
				_back_edges[_edges[edge_idx].get_v2_idx()].push_back(edge_idx);
	}

	// Walk the tree backwards starting at vertex @downstream_v_idx to get a
	// list of edge indices into which the end of the contained read maps.
	edge_idx_t mapped_edges[100];
	size_t num_mapped_edges = 0;
	walk_back_edges(downstream_v_idx, overhang_len,
			mapped_edges, num_mapped_edges, 100);

	// If there was at least one mapped edge but the maximum was not
	// reached, increment the mapped read count on the corresponding edges,
	// but divide by the number of mapped locations to weight the mapping
	// appropriately.
	if (num_mapped_edges < 100 && num_mapped_edges > 0) {
		double div = 1.0f / double(num_mapped_edges);
		for (size_t i = 0; i < num_mapped_edges; i++) {
			assert(mapped_edges[i] < num_edges());
			_edges[mapped_edges[i]].increment_mapped_read_count(div);
		}
	}
}

//
// Given a bidirected string graph, build the corresponding directed string
// graph.
//
void DirectedStringGraph::build_from_bidigraph(const BidirectedStringGraph & bidigraph)
{
	assert(num_vertices() == bidigraph.num_vertices() * 2);

	// Go through each bidirected edge in the bidirected string graph and
	// add the corresponding 2 directed edges to the directed string graph.
	foreach (const BidirectedStringGraphEdge & e, bidigraph.edges()) {
		add_edge_pair(e.get_v1_idx(),
			      e.get_v2_idx(),
			      e.get_dirs(),
			      e.get_seq_1_to_2(), 0, e.length() - 1, false,
			      e.get_seq_2_to_1(), 0, e.length() - 1, false);
	}
	_orig_num_reads = bidigraph._orig_num_reads;
}

void DirectedStringGraph::calculate_A_statistics()
{
	size_t num_reads = _orig_num_reads;
	size_t genome_len;

	size_t bootstrap_genome_len = 0;
	size_t bootstrap_num_reads = num_reads;
	size_t n_edges = num_edges();
	for (size_t i = 0; i < n_edges; i++)
		if (_edges[i].get_v1_idx() & 1)
			bootstrap_genome_len += _edges[i].length();

	float global_arrival_rate = FLOAT_DIV_NONZERO(bootstrap_num_reads,
						      bootstrap_genome_len);

	info("Bootstrapping with bootstream_genome_len = %zu, "
	     "bootstrap_num_reads = %zu, global_arrival_rate = %f",
	     bootstrap_genome_len, bootstrap_num_reads, global_arrival_rate);

	const size_t NUM_BOOTSTRAP_ITERATIONS = 3;
	const float SINGLE_COPY_THRESHOLD = 17.0;

	for (size_t i = 0; i < NUM_BOOTSTRAP_ITERATIONS; i++) {
		size_t num_unique_edges = 0;
		size_t num_optional_edges = 0;
		size_t num_required_edges = 0;

		bootstrap_genome_len = 0;
		bootstrap_num_reads = 0;
		for (size_t j = 0; j < n_edges; j++) {
			DirectedStringGraphEdge & e = _edges[j];
			size_t edge_len = _edges[j].length();
			unsigned edge_reads = _edges[j].get_mapped_read_count();
			float A_statistic = (global_arrival_rate * edge_len) -
					     (edge_reads * M_LN2);
			e.set_A_statistic(A_statistic);
			//info("edge_len = %zu, edge_reads = %u, A_statistic = %f",
			     //edge_len, edge_reads, A_statistic);
			if (A_statistic >= SINGLE_COPY_THRESHOLD) {
				bootstrap_genome_len += edge_len;
				bootstrap_num_reads += edge_reads;
				num_unique_edges++;
			} else if (edge_reads == 0) {
				num_optional_edges++;
			} else {
				num_required_edges++;
			}
		}
		global_arrival_rate = FLOAT_DIV_NONZERO(bootstrap_num_reads,
							bootstrap_genome_len);
		genome_len = DIV_NONZERO(num_reads, global_arrival_rate);
		info("Iteration %zu of %zu:  Estimated genome length "
		     "%zu", i, NUM_BOOTSTRAP_ITERATIONS, genome_len);
		info("num_unique_edges = %zu", num_unique_edges);
		info("num_optional_edges = %zu", num_optional_edges);
		info("num_required_edges = %zu", num_required_edges);
	}
}

void DirectedStringGraph::min_cost_circulation()
{
	info("Initializing lower and upper flow bounds on %zu edges", num_edges());
	foreach (DirectedStringGraphEdge & e, _edges) {
		if (e.get_A_statistic() > 0) {
			e.set_flow_bounds(1, 1);
		} else if (e.get_num_inner_vertices() > 0) {
			e.set_flow_bounds(1, DirectedStringGraphEdge::INFINITE_FLOW);
		} else {
			e.set_flow_bounds(0, DirectedStringGraphEdge::INFINITE_FLOW);
		}
		e.set_cost_per_unit_flow(1);
	}

	v_idx_t n_verts = num_vertices();
	_vertices.resize(n_verts + 2);
	_vertices[n_verts].set_special();
	_vertices[n_verts + 1].set_special();

	info("Adding special vertex and edges");
	for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
		DirectedStringGraphEdge * e;

		//v_idx_t special_v_idx = n_verts | (v_idx & 1);
		for (v_idx_t special_v_idx = n_verts;
		     special_v_idx < n_verts + 2;
		     special_v_idx++)
		{
			e = &add_unlabeled_edge(v_idx, special_v_idx);
			e->set_flow_bounds(0, DirectedStringGraphEdge::INFINITE_FLOW);
			e->set_cost_per_unit_flow(DirectedStringGraphEdge::INFINITE_COST);
			e->set_special();

			e = &add_unlabeled_edge(special_v_idx, v_idx);
			e->set_flow_bounds(0, DirectedStringGraphEdge::INFINITE_FLOW);
			e->set_cost_per_unit_flow(DirectedStringGraphEdge::INFINITE_COST);
			e->set_special();
		}
	}

	n_verts += 2;
	edge_idx_t n_edges = num_edges();

	info("Creating lemon::SmartDigraph with %zu nodes and %zu arcs",
	     n_verts, n_edges);
	lemon::SmartDigraph G;
	G.reserveNode(n_verts);
	G.reserveArc(n_edges);
	lemon::SmartDigraph::ArcMap<int> lower_map(G);
	lemon::SmartDigraph::ArcMap<int> upper_map(G);
	lemon::SmartDigraph::NodeMap<int> supply_map(G);
	for (v_idx_t i = 0; i < n_verts; i++) {
		supply_map[G.addNode()] = 0;
	}
	foreach (DirectedStringGraphEdge & e, _edges) {
		lemon::SmartDigraph::Node node1 = G.nodeFromId(e.get_v1_idx());
		lemon::SmartDigraph::Node node2 = G.nodeFromId(e.get_v2_idx());
		lemon::SmartDigraph::Arc arc = G.addArc(node1, node2);
		lower_map[arc] = e.get_flow_lower_bound();
		upper_map[arc] = e.get_flow_upper_bound();
	}
	info("Initializing network simplex data");

	typedef lemon::NetworkSimplex<lemon::SmartDigraph> simplex_t;
	simplex_t simplex(G);
	simplex.lowerMap(lower_map);
	simplex.upperMap(upper_map);
	simplex.supplyMap(supply_map);
	info("Running network simplex algorithm");
	simplex_t::ProblemType res = simplex.run();
	if (res == simplex_t::INFEASIBLE) {
		fatal_error("No feasible solution to min-cost circulation");
	} else if (res == simplex_t::UNBOUNDED) {
		fatal_error("Objective function unbounded");
	}
	info("Extracting network flow solution");
	for (edge_idx_t edge_idx = 0; edge_idx < n_edges; edge_idx++) {
		lemon::SmartDigraph::Arc arc = G.arcFromId(edge_idx);
		_edges[edge_idx].set_traversal_count(simplex.flow(arc));
	}
	info("Done");
}
