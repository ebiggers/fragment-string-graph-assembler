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
protected:
	float _mapped_read_count;
	unsigned _num_inner_vertices;
	float _A_statistic;
	int _traversal_count;
	bool _is_special;

	StringGraphEdge()
	{
		_mapped_read_count = 1.0;
		_num_inner_vertices = 0;
		_A_statistic = 0.0;
		_traversal_count = 0;
		_is_special = false;
	}

	void print(std::ostream & os) const {
		os << "_mapped_read_count=" << _mapped_read_count << '\t';
		os << "_num_inner_vertices=" << _num_inner_vertices << '\t';
		os << "_A_statistic=" << _A_statistic << '\t';
		os << "_traversal_count=" << _traversal_count << '\t';
		os << "_is_special=" << _is_special << '\t';
	}

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _mapped_read_count;
		ar & _num_inner_vertices;
		ar & _A_statistic;
		ar & _is_special;
		ar & _traversal_count;
	}
public:
	void increment_mapped_read_count(float n) { _mapped_read_count += n; }
	float get_mapped_read_count() const { return _mapped_read_count; }
	void set_mapped_read_count(float n) { _mapped_read_count = n; }

	void increment_num_inner_vertices() { _num_inner_vertices++; }
	unsigned get_num_inner_vertices() const { return _num_inner_vertices; }
	void set_num_inner_vertices(unsigned n) { _num_inner_vertices = n; }

	float get_A_statistic() const { return _A_statistic; }
	void set_A_statistic(float f) { _A_statistic = f; }

	void set_special() { _is_special = true; }
	bool is_special() const { return _is_special; }

	void set_traversal_count(const int traversal_count)
	{
		_traversal_count = traversal_count;
	}
	int get_traversal_count() const { return _traversal_count; }

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
	bool _is_special;

	// Serialize or deserialize the vertex to/from a stream.
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _edge_indices;
		ar & _is_special;
	}

	StringGraphVertex() { _is_special = false; }
public:
	// Add a new edge index to the list of edge indices associated with this
	// string graph vertex.
	void add_edge_idx(const edge_idx_t edge_idx)
	{
		_edge_indices.push_back(edge_idx);
	}

	void set_special() { _is_special = true; }
	bool is_special() const { return _is_special; }

	// Return a const reference to the vertex of edge indices of this string
	// graph vertex.
	const std::vector<edge_idx_t> & edge_indices() const
	{
		return _edge_indices;
	}

	// Return a reference to the vertex of edge indices of this string graph
	// vertex.
	std::vector<edge_idx_t> & edge_indices()
	{
		return _edge_indices;
	}

	edge_idx_t first_edge_idx() const
	{
		assert2(_edge_indices.size() > 0);
		return _edge_indices[0];
	}

	// Print the string graph vertex.
	friend std::ostream & operator<<(std::ostream & os,
					 const StringGraphVertex & v)
	{
		os << "StringGraphVertex {_edge_indices = [";
		foreach(edge_idx_t idx, v.edge_indices())
			os << idx << ", ";
		return os << "] }";
	}
};

template<class VERTEX_t, class EDGE_t, class IMPL_t>
class StringGraph {
protected:
	typedef typename EDGE_t::v_idx_t v_idx_t;
	typedef typename VERTEX_t::edge_idx_t edge_idx_t;

	static const v_idx_t TAG_F_B = 0x0;
	static const v_idx_t TAG_F_E = 0x2;
	static const v_idx_t TAG_G_B = 0x0;
	static const v_idx_t TAG_G_E = 0x1;

	void add_edge_pair(const v_idx_t read_1_idx,
			   const v_idx_t read_2_idx,
			   const v_idx_t dirs,
			   const BaseVec & bv1,
			   const BaseVec::size_type beg_1,
			   const BaseVec::size_type end_1,
			   const bool bv1_rc,
			   const BaseVec & bv2,
			   const BaseVec::size_type beg_2,
			   const BaseVec::size_type end_2,
			   const bool bv2_rc)
	{
		static_cast<IMPL_t*>(this)->add_edge_pair(read_1_idx,
							  read_2_idx,
							  dirs,
							  bv1, beg_1, end_1, bv1_rc,
							  bv2, beg_2, end_2, bv2_rc);
	}
private:

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
		bool rc;

		o.get(f_idx, f_beg, f_end, g_idx, g_beg, g_end, rc);

		const BaseVec & f = bvv[f_idx];
		const BaseVec & g = bvv[g_idx];

		// Skip self-overlaps
		//if (f_idx == g_idx)
			//return;

		// Skip contained overlaps
		if ((f_beg == 0 && f_end == f.size() - 1)
		    || (g_beg == 0 && g_end == g.size() - 1))
			return;

		//info("f_idx = %zu, g_idx = %zu, f_beg=%zu, f_end=%zu",
				//f_idx,g_idx,f_beg,f_end);

		if (f_beg > 0) {
			if (rc) {
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
				 */
				assert2(g_beg > 0);
				assert2(f_beg > 0);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_E | TAG_G_B,
					      g, 0, g_beg - 1, true,
					      f, 0, f_beg - 1, true);
			} else {
				/*
				 *  f.B --------------> f.E
				 *         g.B ----------------> g.E
				 *
				 *  Add f.E -> g.E, g.B -> f.B
				 *
				 *  Or bidirected edge:
				 *
				 *  f >----------> g
				 */

				assert2(g_end + 1 <= g.size() - 1);
				assert2(f_beg > 0);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_E | TAG_G_E,
					      g, g_end + 1, g.size() - 1, false,
					      f, 0, f_beg - 1, true);
			}
		} else {
			if (rc) {
				/*
				 *        f.B ---------------> f.E
				 * g.E <-------------- g.B
				 *
				 *  Add f.B -> g.E, g.B -> f.E
				 *
				 *  Or bidirected edge:
				 *
				 *  f <----------> g
				 */
				assert2(g_end + 1 <= g.size() - 1);
				assert2(f_end + 1 <= f.size() - 1);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_B | TAG_G_E,
					      g, g_end + 1, g.size() - 1, false,
					      f, f_end + 1, f.size() - 1, false);
			} else {

				/*
				 *        f.B ---------------> f.E
				 * g.B --------------> g.E
				 *
				 *  Add f.B -> g.B, g.E -> f.E
				 *
				 *  Or bidirected edge:
				 *
				 *  f <----------< g
				 */

				assert2(g_beg > 0);
				assert2(f_end + 1 <= f.size() - 1);
				add_edge_pair(f_idx, g_idx,
					      TAG_F_B | TAG_G_B,
					      g, 0, g_beg - 1, true,
					      f, f_end + 1, f.size() - 1, false);
			}
		}
	}

protected:
	// Vector of the graph's vertices.
	std::vector<VERTEX_t> _vertices;

	// Vector of the graph's edges.
	std::vector<EDGE_t> _edges;

public:
	size_t _orig_num_reads;
protected:

	// Serialize or deserialize this string graph to/from a stream.
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _vertices;
		ar & _edges;
		ar & _orig_num_reads;
	}

	// Constructor is protected--- use DirectedStringGraph or
	// BidirectedStringGraph instead.
	StringGraph() : _vertices(), _edges()
	{ }

	// Return %true iff the edge structures are large enough to hold indices
	// for @num_vertices_needed different vertices.
	bool enough_v_indices(size_t num_vertices_needed) const
	{
		return num_vertices_needed <= std::numeric_limits<v_idx_t>::max();
	}

	// Add an edge to the vector of edges of this string graph and return
	// its index.
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
		foreach(VERTEX_t & v, _vertices)
			std::sort(v.edge_indices().begin(), v.edge_indices().end(), cmp);
	}
public:

	// Return a reference to a vector of this string graph's edges.
	std::vector<EDGE_t> & edges() { return _edges; }
	const std::vector<EDGE_t> & edges() const { return _edges; }

	// Return a reference to a vector of this string graph's vertices.
	std::vector<VERTEX_t> & vertices() { return _vertices; }
	const std::vector<VERTEX_t> & vertices() const { return _vertices; }

	// Return the number of edges in this string graph.
	edge_idx_t num_edges() const { return _edges.size(); }

	// Return the number of vertices in this string graph.
	v_idx_t num_vertices() const { return _vertices.size(); }

	// Return the index of the edge f -> g, which must exist */
	edge_idx_t locate_edge(const v_idx_t f_idx, const v_idx_t g_idx) const
	{
		assert(f_idx < num_vertices() && g_idx < num_vertices());
		foreach(edge_idx_t edge_idx, _vertices[f_idx].edge_indices()) {
			assert(edge_idx < num_edges());
			if (_edges[edge_idx].get_v2_idx() == g_idx)
				return edge_idx;
		}
		assert(0);
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
		char buf[10];
		in.read(buf, 10);
		const char * magic = static_cast<IMPL_t*>(this)->magic;
		if (memcmp(buf, magic, 10) != 0) {
			// Need to throw exception because of AnyStringGraph.h
			throw std::runtime_error("Invalid magic characters in graph file");
		}
		boost::archive::binary_iarchive ar(in);
		ar >> *this;
		static_cast<const IMPL_t*>(this)->assert_graph_valid();
	}

	// Write this string graph to a file.
	void write(const char *filename) const
	{
		static_cast<const IMPL_t*>(this)->assert_graph_valid();
		std::ofstream out(filename);
		out.write(static_cast<const IMPL_t*>(this)->magic, 10);
		boost::archive::binary_oarchive ar(out);
		ar << *this;
		out.close();
		if (!out)
			fatal_error_with_errno("Error writing to \"%s\"", filename);
	}

	// Print this string graph.
	void print(std::ostream & os, const bool print_seqs) const
	{
		for (v_idx_t v_idx = 0; v_idx < _vertices.size(); v_idx++) {
			const VERTEX_t & v = _vertices[v_idx];
			foreach(const edge_idx_t edge_idx, v.edge_indices()) {
				_edges[edge_idx].print(os, v_idx, print_seqs);
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
	void print_dot(std::ostream & os, const bool print_seqs) const
	{
		os << "digraph {\n"
		   << "\tnode [shape=circle fontname=\"Arial\"]\n"
		   << "\tedge [fontname=\"Courier new bold\" fontsize=11]\n";
		print_dot_graph_attribs(os);
		for (v_idx_t v_idx = 0; v_idx < _vertices.size(); v_idx++)
			_vertices[v_idx].print_dot(os, v_idx);

		for (v_idx_t v_idx = 0; v_idx < _vertices.size(); v_idx++)
			foreach(edge_idx_t edge_idx, _vertices[v_idx].edge_indices())
				_edges[edge_idx].print_dot(os, v_idx, print_seqs);
		os << "}" << std::endl;
	}

	// Builds this string graph from a set of reads and their overlaps.
	void build(const BaseVecVec & bvv, const OverlapVecVec & ovv)
	{
		assert(bvv.size() == ovv.size());
		foreach(const OverlapVecVec::OverlapSet & overlap_set, ovv) {
			foreach(const Overlap & o, overlap_set) {
				assert_overlap_valid(o, bvv, 1, 0);
				add_edge_from_overlap(bvv, o);
			}
		}
		info("String graph has %zu vertices and %zu edges",
		     num_vertices(), num_edges());
		info("Average of %.2f edges per vertex",
		     DOUBLE_DIV_NONZERO(num_edges(), num_vertices()));
	}
};
