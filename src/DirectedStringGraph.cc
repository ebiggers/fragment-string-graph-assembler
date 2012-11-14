#include "DirectedStringGraph.h"
#include <fstream>
#include <ostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

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
		// map from the adjacent vertices' indices to the back edges
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
						 const std::vector<bool> & v_inner)
{
	BaseVecVec::size_type new_seq_len = e.length();
	v_idx_t vi_idx = e.get_v2_idx();
	assert(v_inner[vi_idx]);

	//size_t path_len = 1;
	//info("v1 = %zu", e.get_v1_idx());
	//info("v2 = %zu", e.get_v2_idx());
	// Found beginning of unbranched path.  Walk along it until
	// the end to get the total sequence length.
	do {
		const DirectedStringGraphVertex &vi = _vertices[vi_idx];
		assert(vi.out_degree() == 1);
		//path_len++;
		const DirectedStringGraphEdge &ei_i1 = _edges[vi.first_edge_idx()];
		if (new_seq_len + ei_i1.length() < new_seq_len)
			fatal_error("Edge too long");
		new_seq_len += ei_i1.length();
		vi_idx = ei_i1.get_v2_idx();
		//info("vi = %zu", vi_idx);
	} while (v_inner[vi_idx]);

	//info("path_len = %zu", path_len);

	BaseVec & new_seq = e.get_seq();
	BaseVec::size_type seq_idx = new_seq.length();
	new_seq.resize(new_seq_len);
	vi_idx = e.get_v2_idx();
	do {
		const DirectedStringGraphVertex &vi = _vertices[vi_idx];
		const edge_idx_t ei_i1_idx = vi.first_edge_idx();
		const DirectedStringGraphEdge &ei_i1 = _edges[ei_i1_idx];
		const BaseVec & ei_i1_seq = ei_i1.get_seq();
		for (BaseVec::size_type i = 0; i < ei_i1_seq.length(); i++) {
			assert2(seq_idx < new_seq_len);
			new_seq.set(seq_idx++, ei_i1_seq[i]);
		}
		e.increment_mapped_read_count(ei_i1.get_mapped_count());
		remove_edge[ei_i1_idx] = true;
		vi_idx = ei_i1.get_v2_idx();
	} while (v_inner[vi_idx]);
	assert(new_seq_len == seq_idx);
	e.set_v2_idx(vi_idx);
}

void DirectedStringGraph::collapse_unbranched_paths()
{
	const v_idx_t n_verts = num_vertices();
	const edge_idx_t n_edges = num_edges();

	info("Collapsing unbranched paths in directed string graph");
	info("Original graph has %u vertices and %u edges", n_verts, n_edges);

	size_t num_inner_vertices = 0;
	// Find whether each vertex is inner or not.  A vertex is inner iff it
	// has indegree 1 and outdegree 1.
	std::vector<bool> v_inner(n_verts, false);
	{
		std::vector<unsigned char> v_in_degrees(n_verts, 0);
		std::vector<unsigned char> v_out_degrees(n_verts, 0);
		for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
			const DirectedStringGraphVertex & v = _vertices[v_idx];
			foreach(edge_idx_t edge_idx, v.edge_indices()) {
				const DirectedStringGraphEdge & e = _edges[edge_idx];
				v_idx_t v1_idx, v2_idx;
				e.get_v_indices(v1_idx, v2_idx);
				if (v_out_degrees[v1_idx] < 2)
					v_out_degrees[v1_idx]++;
				if (v_in_degrees[v2_idx] < 2)
					v_in_degrees[v2_idx]++;
			}
		}
		for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
			if (v_in_degrees[v_idx] == 1 && v_out_degrees[v_idx] == 1) {
				v_inner[v_idx] = true;
				num_inner_vertices++;
			}
		}
	}

	info("Found %zu inner vertices (%f%% of all vertices)",
	     num_inner_vertices,
	     (n_verts ? 100 * double(num_inner_vertices) / n_verts : 0));

	// Go through each non-inner vertex and look for any neighboring inner
	// vertices.  These are the starts of unbranched paths that will be
	// collapsed.  For each such path, follow it until its end and determine
	// the total length of the edges.  Then, set the first edge's label to
	// the concatenation of the labels of all the edges in the unbranched
	// path, then make the first edge point directly to the last vertex in
	// the path and mark the other edges for deletion.
	size_t num_unbranched_paths = 0;
	std::vector<bool> remove_edge(n_edges, false);
	for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
		if (!v_inner[v_idx]) {
			const DirectedStringGraphVertex & v = _vertices[v_idx];
			foreach(edge_idx_t edge_idx, v.edge_indices()) {
				DirectedStringGraphEdge & e = _edges[edge_idx];
				if (v_inner[e.get_v2_idx()]) {
					num_unbranched_paths++;
					follow_unbranched_path(e, remove_edge, v_inner);
				}
			}
		}
	}

	info("Found %zu unbranched paths", num_unbranched_paths);

	// Compute the new vertex indices.
	std::vector<v_idx_t> new_v_indices(n_verts);
	v_idx_t new_v_idx = 0;
	for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++)
		if (v_inner[v_idx])
			new_v_indices[v_idx] = std::numeric_limits<v_idx_t>::max();
		else
			new_v_indices[v_idx] = new_v_idx++;

	info("Updated vertices are indexed [0, %lu)", new_v_idx);

	info("Updating edges");
	// Compute the new edge indices, set the new vertex indices in each edge,
	// and move the edges
	std::vector<edge_idx_t> new_edge_indices(n_edges);
	edge_idx_t new_edge_idx = 0;
	for (edge_idx_t edge_idx = 0; edge_idx < n_edges; edge_idx++) {
		if (!remove_edge[edge_idx]) {
			DirectedStringGraphEdge & e = _edges[edge_idx];
			new_edge_indices[edge_idx] = new_edge_idx;
			v_idx_t v1_idx, v2_idx;
			e.get_v_indices(v1_idx, v2_idx);
			e.set_v_indices(new_v_indices[v1_idx], new_v_indices[v2_idx]);
			_edges[new_edge_idx++] = _edges[edge_idx];
		}
	}
	info("Updated edges are indexed [0, %lu)", new_edge_idx);
	info("%zu edges were removed (%f%% of total)",
	     _edges.size() - new_edge_idx,
	     (_edges.size() ? 100.0 * (_edges.size() - new_edge_idx) / _edges.size() : 0));

	_edges.resize(new_edge_idx);

	info("Updating vertices");
	// Set new edge indices in each vertex and move the vertices
	new_v_idx = 0;
	for (v_idx_t v_idx = 0; v_idx < n_verts; v_idx++) {
		if (!v_inner[v_idx]) {
			DirectedStringGraphVertex & v = _vertices[v_idx];
			std::vector<edge_idx_t> & edge_indices = v.edge_indices();
			for (size_t i = 0; i < v.out_degree(); i++) {
				edge_idx_t edge_idx = edge_indices[i];
				edge_idx_t new_edge_idx = new_edge_indices[edge_idx];
				edge_indices[i] = new_edge_idx;
			}
			_vertices[new_v_idx++] = _vertices[v_idx];
		}
	}
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
		assert(e.get_mapped_count() >= 1.0);
		total_mapped_count += e.get_mapped_count();
		if (e.get_mapped_count() > 1.0)
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

	v_idx_t num_components = 0;
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
	os << "    Number of components: " << num_components << std::endl;
	os << "    Component sizes:" << std::endl;
	foreach (v_idx_t component_size, component_sizes) {
		os << "        " << component_size << ' '
		   << (component_size == 1 ? "vertex" : "vertices") << std::endl;
	}
	os << "}" << std::endl;
}

void DirectedStringGraph::walk_back_edges(const v_idx_t v_idx,
					  size_t overhang_len,
					  edge_idx_t mapped_edges[],
					  size_t & num_mapped_edges,
					  size_t max_mapped_edges)
{
	assert(_back_edges.size() == _vertices.size());

	foreach (edge_idx_t edge_idx, _back_edges[v_idx]) {
		const DirectedStringGraphEdge & e = _edges[edge_idx];
		assert(e.get_v2_idx() == v_idx);
		if (overhang_len < e.length()) {
			mapped_edges[num_mapped_edges++] = edge_idx;
			if (num_mapped_edges == max_mapped_edges)
				return;
		}
	}

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

void DirectedStringGraph::map_contained_read(size_t contained_read_idx,
					     const Overlap & o,
					     const size_t overhang_len)
{
	Overlap::read_idx_t f_idx;
	Overlap::read_pos_t f_beg;
	Overlap::read_pos_t f_end;
	Overlap::read_idx_t g_idx;
	Overlap::read_pos_t g_beg;
	Overlap::read_pos_t g_end;
	bool rc;

	o.get(f_idx, f_beg, f_end, g_idx, g_beg, g_end, rc);

	assert(contained_read_idx < num_vertices());

	if (contained_read_idx == g_idx) {
		std::swap(f_idx, g_idx);
		std::swap(f_beg, g_beg);
		std::swap(f_end, g_end);
	} else {
		assert(contained_read_idx == f_idx);
	}

	assert(f_idx * 2 + 1 < num_vertices());
	assert(g_idx * 2 + 1 < num_vertices());

	v_idx_t v_idx = contained_read_idx * 2;

	if (rc) {
		assert(overhang_len == g_beg);
	} else {
		// overhang_len == g.length() - 1 - g_end;
		//
		// but we don't have g here.
		v_idx++;
	}

	if (_back_edges.size() == 0) {
		info("Indexing back edges (num_vertices = %zu, num_edges = %zu)",
				num_vertices(), num_edges());
		_back_edges.resize(num_vertices());
		foreach (const DirectedStringGraphVertex v, _vertices)
			foreach (edge_idx_t edge_idx, v.edge_indices())
				_back_edges[_edges[edge_idx].get_v2_idx()].push_back(edge_idx);
	}
	assert(_back_edges.size() == num_vertices());

	assert(v_idx < num_vertices());

	edge_idx_t mapped_edges[100];
	size_t num_mapped_edges = 0;
	walk_back_edges(v_idx, overhang_len, mapped_edges, num_mapped_edges, 100);

	if (num_mapped_edges < 100 && num_mapped_edges > 0) {
		double div = 1.0f / double(num_mapped_edges);
		for (size_t i = 0; i < num_mapped_edges; i++) {
			assert(mapped_edges[i] < num_edges());
			_edges[mapped_edges[i]].increment_mapped_read_count(div);
		}
	}
}
