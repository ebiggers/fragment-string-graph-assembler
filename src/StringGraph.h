#pragma once

#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "util.h"
#include <fstream>

#include "Overlap.h"
#include "BaseVecVec.h"

// Base class for edges of the string graph.
class StringGraphEdge {
};

// Base class for vertices of the string graph.
class StringGraphVertex {
public:
	// Unsigned integer type of an edge index.  This places one upper bound
	// an the number of edges that can be in the graph.
	typedef unsigned long edge_idx_t;
protected:
	// vector of the indices of the edges that go out from this vertex.
	std::vector<edge_idx_t> _edge_indices;

	// Serialize or deserialize the vertex to/from a stream.
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _edge_indices;
	}

	StringGraphVertex() { }
public:
	// Add a new edge index to the list of edge indices associated with this
	// string graph vertex.
	void add_edge_idx(const edge_idx_t edge_idx)
	{
		_edge_indices.push_back(edge_idx);
	}

	// Return a const reference to the vertex of edge indices of this string
	// graph vertex.
	const std::vector<edge_idx_t> & edge_indices() const
	{
		return _edge_indices;
	}

	edge_idx_t first_edge_idx() const
	{
		return _edge_indices[0];
	}

	// Return a reference to the vertex of edge indices of this string graph
	// vertex.
	std::vector<edge_idx_t> & edge_indices()
	{
		return _edge_indices;
	}

	// Print the string graph vertex.
	friend std::ostream & operator<<(std::ostream & os,
					 const StringGraphVertex & v)
	{
		os << "StringGraphVertex {_edge_indices = [";
		for (edge_idx_t idx : v.edge_indices())
			os << idx << ", ";
		return os << "] }";
	}
};

template<class VERTEX_t, class EDGE_t, class IMPL_t>
class StringGraph {
protected:
	typedef typename EDGE_t::v_idx_t v_idx_t;
	typedef typename VERTEX_t::edge_idx_t edge_idx_t;

	void add_edge_pair(const v_idx_t read_1_idx,
			   const v_idx_t read_2_idx,
			   const v_idx_t dirs,
			   const BaseVec & bv1,
			   const BaseVec::size_type beg_1,
			   const BaseVec::size_type end_1,
			   const BaseVec & bv2,
			   const BaseVec::size_type beg_2,
			   const BaseVec::size_type end_2)
	{
		static_cast<IMPL_t*>(this)->add_edge_pair(read_1_idx,
							  read_2_idx,
							  dirs,
							  bv1, beg_1, end_1,
							  bv2, beg_2, end_2);
	}
private:

	static const v_idx_t TAG_F_B = 0x0;
	static const v_idx_t TAG_F_E = 0x2;
	static const v_idx_t TAG_G_B = 0x0;
	static const v_idx_t TAG_G_E = 0x1;

	// Given an uncontained overlap and the read set from which it came, add
	// the corresponding edge(s) to this string graph.
	void add_edge_from_overlap(const BaseVecVec & bvv, const Overlap & o)
	{
		Overlap::read_idx_t f_idx;
		Overlap::read_pos_t f_beg;
		Overlap::read_pos_t f_end;
		Overlap::read_idx_t g_idx;
		Overlap::read_pos_t g_beg;
		Overlap::read_pos_t g_end;

		o.get(f_idx, f_beg, f_end, g_idx, g_beg, g_end);

		const BaseVec & f = bvv[f_idx];
		const BaseVec & g = bvv[g_idx];

		// Skip self-overlaps
		//if (f_idx == g_idx)
			//return;

		// Skip contained overlaps
		if (((f_beg == 0 && f_end == f.size() - 1)
		    || (f_beg == f.size() - 1 && f_end == 0)
		    || (g_beg == 0 && g_end == g.size() - 1)
		    || (g_beg == g.size() - 1 && g_end == 0)))
			return;

		if (f_beg > f_end) {
			std::swap(f_idx, g_idx);
			std::swap(f_beg, g_beg);
			std::swap(f_end, g_end);
		}
		if (f_beg > 0) {
			if (g_beg < g_end) {
				/*
				 *  f.B --------------> f.E
				 *         g.B ----------------> g.E
				 *
				 *  Add f.E -> g.E, g.B -> f.B
				 *
				 *  Or bidirected edge:
				 *
				 *  f >----------> g
				 *
				 *     f.E -> g.E label: g[g_end + 1 ... g.size() - 1]
				 *     g.B -> f.B label: f[0 ... f_beg - 1]
				 */

				assert2(g_end + 1 <= g.size() - 1);
				assert2(f_beg > 0);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_E | TAG_G_E,
					      g, g_end + 1, g.size() - 1,
					      f, f_beg - 1, 0);
			} else {

				/*
				 *  f.B --------------> f.E
				 *         g.E <---------------  g.B
				 *
				 *  Add f.E -> g.B, g.E -> f.B
				 *
				 *  Or bidirected edge:
				 *
				 *  f >----------< g
				 *
				 *     f.E -> g.B label: g[g_end - 1 ... 0]
				 *     g.E -> f.B label: f[f_beg - 1 ... 0]
				 */

				assert2(g_end > 0);
				assert2(f_beg > 0);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_E | TAG_G_B,
					      g, g_end - 1, 0,
					      f, f_beg - 1, 0);
			}
		} else {
			if (g_beg < g_end) {

				/*
				 *        f.B ---------------> f.E
				 * g.B --------------> g.E
				 *
				 *  Add f.B -> g.B, g.E -> f.E
				 *
				 *  Or bidirected edge:
				 *
				 *  f <----------< g
				 *
				 *     f.B -> g.B label: g[g_beg - 1 ... 0]
				 *     g.E -> f.E label: f[f_end + 1 ... f.size() - 1]
				 */

				assert2(g_beg > 0);
				assert2(f_end + 1 <= f.size() - 1);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_B | TAG_G_B,
					      g, g_beg - 1, 0,
					      f, f_end + 1, f.size() - 1);
			} else {

				/*
				 *        f.B ---------------> f.E
				 * g.E <-------------- g.B
				 *
				 *  Add f.B -> g.E, g.B -> f.E
				 *
				 *  Or bidirected edge:
				 *
				 *  f <----------> g
				 *
				 *     f.B -> g.E label: g[g_beg + 1 ... g.size() - 1]
				 *     g.B -> f.E label: f[f_end + 1 ... f.size() - 1]
				 */
				assert2(g_beg + 1 <= g.size() - 1);
				assert2(f_end + 1 <= f.size() - 1);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_B | TAG_G_E,
					      g, g_beg + 1, g.size() - 1,
					      f, f_end + 1, f.size() - 1);
			}
		}
	}

protected:
	// Vector of the graph's vertices.
	std::vector<VERTEX_t> _vertices;

	// Vector of the graph's edges.
	std::vector<EDGE_t> _edges;

	// Serialize or deserialize this string graph to/from a stream.
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _vertices;
		ar & _edges;
	}

	// Constructor is protected--- use DirectedStringGraph or
	// BidirectedStringGraph instead.
	StringGraph() { }

	// Return %true iff the edge structures are large enough to hold indices
	// for @num_vertices_needed different vertices.
	bool enough_v_indices(size_t num_vertices_needed) const
	{
		return num_vertices_needed <= std::numeric_limits<v_idx_t>::max();
	}

	// Add an edge to vector of edges of this string graph and return its
	// index.
	edge_idx_t push_back_edge(const EDGE_t & e)
	{
		if (_edges.size() >= std::numeric_limits<edge_idx_t>::max())
			fatal_error("Too many edges");
		edge_idx_t edge_idx = _edges.size();
		_edges.push_back(e);
		return edge_idx;
	}

private:
	class cmp_by_edge_length {
	private:
		const std::vector<EDGE_t> & _edges;
	public:
		cmp_by_edge_length(const std::vector<EDGE_t> & edges)
			: _edges(edges) { }

		template <typename edge_idx_t>
		bool operator()(edge_idx_t edge_idx_1, edge_idx_t edge_idx_2) const
		{
			return _edges[edge_idx_1].length() < _edges[edge_idx_2].length();
		}
	};
protected:
	void sort_adjlists_by_edge_len()
	{
		cmp_by_edge_length cmp(_edges);
		for (VERTEX_t & v : _vertices)
			std::sort(v.edge_indices().begin(), v.edge_indices().end(), cmp);
	}
public:

	// Return a reference to a vector of this string graph's edges.
	std::vector<EDGE_t> & edges()
	{
		return _edges;
	}

	// Return a reference to a vector of this string graph's vertices.
	std::vector<VERTEX_t> & vertices()
	{
		return _vertices;
	}

	// Return the number of edges in this string graph.
	size_t num_edges() const
	{
		return _edges.size();
	}

	// Return the number of vertices in this string graph.
	size_t num_vertices() const
	{
		return _vertices.size();
	}

	// Delete all edges and vertices from this string graph.
	void clear()
	{
		_edges.resize(0);
		_vertices.resize(0);
	}

	// Read this string graph from a file.
	void read(const char *filename)
	{
		this->clear();
		std::ifstream in(filename);
		if (!in)
			fatal_error_with_errno("Error opening \"%s\"", filename);
		boost::archive::binary_iarchive ar(in);
		ar >> *this;
	}

	// Write this string graph to a file.
	void write(const char *filename) const
	{
		std::ofstream out(filename);
		boost::archive::binary_oarchive ar(out);
		ar << *this;
		out.close();
		if (out.bad())
			fatal_error_with_errno("Error writing to \"%s\"", filename);
	}

	// Print this string graph.
	void print(std::ostream & os) const
	{
		for (v_idx_t v_idx = 0; v_idx < _vertices.size(); v_idx++) {
			const VERTEX_t & v = _vertices[v_idx];
			for (const edge_idx_t edge_idx : v.edge_indices()) {
				_edges[edge_idx].print(os, v_idx);
				os << '\n';
			}
		}
		os << std::flush;
	}

	void print_dot_graph_attribs(std::ostream & os) const
	{
		static_cast<const IMPL_t*>(this)->print_dot_graph_attribs(os);
	}

	// Print this string graph in DOT format.
	void print_dot(std::ostream & os) const
	{
		os << "digraph {\n";
		print_dot_graph_attribs(os);
		for (v_idx_t v_idx = 0; v_idx < _vertices.size(); v_idx++)
			_vertices[v_idx].print_dot(os, v_idx);

		for (v_idx_t v_idx = 0; v_idx < _vertices.size(); v_idx++)
			for (edge_idx_t edge_idx : _vertices[v_idx].edge_indices())
				_edges[edge_idx].print_dot(os, v_idx);
		os << "}" << std::endl;
	}

	// Builds this string graph from a set of reads and their overlaps.
	void build(const BaseVecVec & bvv, const OverlapVecVec & ovv)
	{
		assert(bvv.size() == ovv.size());
		for (auto overlap_set : ovv) {
			for (const Overlap & o : overlap_set) {
				assert_overlap_valid(o, bvv, 0, 0);
				add_edge_from_overlap(bvv, o);
			}
		}
		info("String graph has %zu vertices and %zu edges",
		     num_vertices(), num_edges());
		info("Average of %.2f edges per vertex",
		     num_vertices() ?
			double(num_edges()) / double(num_vertices()) : 0);
	}
};
