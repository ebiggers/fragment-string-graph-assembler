#include "DirectedStringGraph.h"
#include <fstream>
#include <ostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

//
// Perform transitive edge reduction on a directed string graph:
// Remove all edges v -> x where there exist edges v -> w -> x.
//
// See:
//    "The fragment string assembly graph", Eugene W. Myers, 2005
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

	// Per-edge boolean indicating whether each edge is to be deleted as
	// part of the transitive reduction or not.
	std::vector<bool> reduce_edge(this->num_edges(), false);

	std::vector<edge_idx_t> in_edges;
	std::vector<unsigned char> num_in_edges;
	std::vector<uint16_t> v_idx_to_in_edges(this->num_vertices(), 0xffff);

	// Iterate through every vertex @v in the graph that has outgoing edges.
	for (size_t v_idx = 0; v_idx < vertices.size(); v_idx++) {
		const DirectedStringGraphVertex & v = vertices[v_idx];

		if (v.out_degree() == 0)
			continue;

		if (v.out_degree() >= 0xffff) {
			fatal_error("Out degree of vertex %zu is too high (%zu, need < 65535)",
				    v_idx, v.out_degree());
		}

		in_edges.clear();
		num_in_edges.clear();
		std::fill(num_in_edges.begin(), num_in_edges.end(), 0);
		uint16_t cur_in_edges_idx = 0;

		// Make it possible to look up edges v -> w given a vertex w.
		// (v is vertices[v_idx]).
		// It's assumed there can be at most 2 edges v -> w, which is
		// true if the edges originate from non-redundant overlaps.
		for (size_t i = 0; i < v.out_degree(); i++) {
			const edge_idx_t edge_idx = v.edge_indices()[i];
			const v_idx_t v2_idx = edges[edge_idx].get_v2_idx();
			uint16_t in_edges_idx = v_idx_to_in_edges[v2_idx];

			//vertex_marks[v2_idx] = INPLAY;

			if (in_edges_idx == 0xffff) {
				in_edges_idx = cur_in_edges_idx++;
				v_idx_to_in_edges[v2_idx] = in_edges_idx;
				in_edges.push_back(edge_idx);
				in_edges.push_back(0);
				num_in_edges.push_back(1);
			} else {
				assert(num_in_edges[in_edges_idx] == 1);
				in_edges[in_edges_idx * 2 + 1] = edge_idx;
				num_in_edges[in_edges_idx] = 2;
			}
		}

		// Length of the longest sequence label on the edges leaving
		// vertex @v.
		const BaseVec::size_type longest = edges[v.edge_indices().back()].length();

		// For each outgoing edge from v -> w in order of labeled
		// sequence length, consider each vertex w that is still marked
		// INPLAY.
		for (const edge_idx_t v_edge_idx : v.edge_indices()) {
			const DirectedStringGraphEdge & edge_v_w = edges[v_edge_idx];
			const v_idx_t w_idx = edge_v_w.get_v2_idx();
			const DirectedStringGraphVertex & w = vertices[w_idx];

			for (const edge_idx_t w_edge_idx : w.edge_indices()) {
				const DirectedStringGraphEdge & edge_w_x = edges[w_edge_idx];

				BaseVecVec::size_type len = edge_v_w.length() +
							    edge_w_x.length();
				if (len > longest)
					continue;

				const v_idx_t x_idx = edge_w_x.get_v2_idx();

				uint16_t x_in_edges_idx = v_idx_to_in_edges[x_idx];

				for (unsigned char in_edge_i = 0;
				     in_edge_i < num_in_edges[x_in_edges_idx];
				     in_edge_i++)
				{
					const edge_idx_t in_edge_idx =
						in_edges[x_in_edges_idx * 2 + in_edge_i];
					const DirectedStringGraphEdge & in_edge = 
						edges[in_edge_idx];

					if (in_edge.length() == 

				}

				if (vertex_marks[x_idx] == INPLAY)
					vertex_marks[x_idx] = ELIMINATED;
			}
		}

		// Once again, go through the outgoing edges from v.  For each
		// neighboring vertex marked ELIMINATED, mark the corresponding
		// edge(s) for reduction.  Return both INPLAY and ELIMINATED
		// vertices to VACANT status.
		for (const edge_idx_t edge_idx : v.edge_indices())
			if (vertex_marks[edges[edge_idx].get_v2_idx()] == ELIMINATED)
				reduce_edge[edge_idx] = true;

		for (const edge_idx_t edge_idx : v.edge_indices())
			vertex_marks[edges[edge_idx].get_v2_idx()] = VACANT;
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
	for (DirectedStringGraphVertex & v : vertices) {
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
			for (edge_idx_t edge_idx : v.edge_indices()) {
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
			for (edge_idx_t edge_idx : v.edge_indices()) {
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
