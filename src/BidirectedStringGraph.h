#pragma once

#include "StringGraph.h"
#include "BaseVec.h"
#include <boost/serialization/access.hpp>
#include <ostream>
#include <inttypes.h>

// An edge of a bidirected string graph.
//
// An edge in a bidirected graph has a head at each edge of the edge.  As a
// result, there are three different orientations:
//
//    1 >----------> 2    (kind of like a normal directed edge)
//    1 >----------< 2    (inward)
//    1 <----------> 2    (outward)
//
// (The fourth possibility is actually just the first flipped around.)
//    1 <----------< 2
//
// Also, since this is a *string* graph, every edge is labeled with a string
// (DNA sequence in this case).  And since this is a bidirected graph, there is
// actually a different string for each way the edge can be traversed.
//
class BidirectedStringGraphEdge : public StringGraphEdge {
private:
	// (high to low bit)
	// 1 bit:   1 iff the head at vertex 1 is directed away from it.
	// 1 bit:   1 iff the head at vertex 2 is directed towards it.
	// 31 bits: index of vertex 1.
	// 31 bits: index of vertex 2.
	uint64_t _data;

	// Sequence when this bidirected edge is traversed in the direction
	// vertex 1 to vertex 2.
	BaseVec _seq_1_to_2;

	// Sequence when this bidirected edge is traversed in the direction
	// vertex 2 to vertex 1.
	BaseVec _seq_2_to_1;

	// Serialize or deserialize this bidirected string graph edge to a
	// stream.
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _data;
		ar & _seq_1_to_2;
		ar & _seq_2_to_1;
	}
public:
	typedef unsigned int v_idx_t;

	BaseVec & get_seq_1_to_2()
	{
		return _seq_1_to_2;
	}
	BaseVec & get_seq_2_to_1()
	{
		return _seq_2_to_1;
	}

	const BaseVec & get_seq_1_to_2() const
	{
		return _seq_1_to_2;
	}

	const BaseVec & get_seq_2_to_1() const
	{
		return _seq_2_to_1;
	}

	BaseVec::size_type length() const
	{
		return _seq_1_to_2.size();
	}

	v_idx_t get_other_v_idx(v_idx_t this_v_idx) const
	{
		v_idx_t v1_idx, v2_idx;
		get_v_indices(v1_idx, v2_idx);
		if (this_v_idx == v1_idx)
			return v2_idx;
		else {
			assert(this_v_idx == v2_idx);
			return v1_idx;
		}
	}

	v_idx_t get_v1_idx() const
	{
		return (_data & 0x7fffffffULL);
	}

	v_idx_t get_v2_idx() const
	{
		return (_data >> 31) & 0x7fffffffULL;
	}

	void set_v1_idx(v_idx_t v1_idx)
	{
		_data &= ~0x7fffffffULL;
		_data |= v1_idx;
	}

	void set_v2_idx(v_idx_t v2_idx)
	{
		_data &= ~(0x7fffffffULL << 31);
		_data |= (uint64_t(v2_idx) << 31);
	}

	void set_v_indices(v_idx_t v1_idx, v_idx_t v2_idx)
	{
		_data = (_data & (3ULL << 62)) | (uint64_t(v2_idx) << 31) | v1_idx;
	}

	void get_v_indices(v_idx_t & v1_idx, v_idx_t & v2_idx) const
	{
		v1_idx = get_v1_idx();
		v2_idx = get_v2_idx();
	}

	v_idx_t get_dirs() const
	{
		return _data >> 62;
	}

	bool v1_outward() const { return (get_dirs() & 0x2) != 0; }
	bool v2_inward() const { return (get_dirs() & 0x1) != 0; }
	bool v1_inward() const { return !(v1_outward()); }
	bool v2_outward() const { return !(v2_inward()); }

	bool this_v_outward(v_idx_t this_v_idx) const
	{
		v_idx_t v1_idx, v2_idx;
		get_v_indices(v1_idx, v2_idx);
		if (this_v_idx == v1_idx) {
			return v1_outward();
		} else {
			assert(this_v_idx == v2_idx);
			return v2_outward();
		}
	}

	bool other_v_outward(v_idx_t this_v_idx) const
	{
		v_idx_t v1_idx, v2_idx;
		get_v_indices(v1_idx, v2_idx);
		if (this_v_idx == v1_idx) {
			return v2_outward();
		} else {
			assert(this_v_idx == v2_idx);
			return v1_outward();
		}
	}


	void set_dirs(unsigned dirs)
	{
		_data = (_data & ~(3ULL << 62)) | (uint64_t(dirs) << 62);
	}

	// Print a bidirected string graph edge
	void print(std::ostream & os, const v_idx_t v_idx) const
	{
		v_idx_t read_1_idx = get_v1_idx();
		v_idx_t read_2_idx = get_v2_idx();
		char head_1 = v1_outward() ? '>' : '<';
		char head_2 = v2_inward() ? '>' : '<';
		const BaseVec * seq = &get_seq_1_to_2();
		if (read_1_idx != v_idx) {
			std::swap(read_1_idx, read_2_idx);
			head_1 = (head_1 == '>') ? '<' : '>';
			head_2 = (head_2 == '>') ? '<' : '>';
			seq = &get_seq_2_to_1();
		}
		os << (read_1_idx + 1) << ' ' << head_1 << "---------"
		   << head_2 << ' ' << (read_2_idx + 1)
		   << "\t\"" << get_seq_1_to_2() << '"';
	}

	// Print a bidirected string graph edge in DOT format
	void print_dot(std::ostream & os, const v_idx_t v_idx) const
	{
		if (v_idx == get_v1_idx()) {
			const char *head_1 = (v1_inward()) ? "normal" : "inv";
			const char *head_2 = (v2_inward()) ? "normal" : "inv";
			os << "\tv" << get_v1_idx()
			   << " -> "
			   << "v" << get_v2_idx()
			   << " ["
			   << " dir=both arrowhead=" << head_2
			   << " arrowtail=" << head_1
			   << " label=\"" << length() << "\""
			   << " ];\n";
		}
	}
};


// A vertex of a bidirected string graph.
class BidirectedStringGraphVertex : public StringGraphVertex {
public:
	size_t out_degree() const
	{
		unimplemented();
	}

	// Print a bidirected string graph vertex in DOT format
	void print_dot(std::ostream & os, size_t v_idx) const
	{
		os << "\tv" << v_idx
		   << " [ label = \"" << (v_idx + 1) << "\" ];\n";
	}
};

// A bidirected string graph.
class BidirectedStringGraph : public StringGraph<BidirectedStringGraphVertex,
						 BidirectedStringGraphEdge,
						 BidirectedStringGraph>
{
public:
	// Initialize a bidirected string graph with enough space for @num_reads
	// reads to be inserted
	BidirectedStringGraph(size_t num_reads)
	{
		if (!enough_v_indices(num_reads * 2))
			fatal_error("Too many reads (%zu)", num_reads);
		_vertices.resize(num_reads);
	}

	// Read a bidirected string graph from a file
	BidirectedStringGraph(const char *filename)
	{
		this->read(filename);
	}


	void print_dot_graph_attribs(std::ostream & os) const
	{
		//os << "\tconcentrate=true;\n";
		//os << "\tnode [shape = rect];\n";
	}

	void transitive_reduction();
	void collapse_unbranched_paths();

	// An an edge to the bidirected string graph, produced from an overlap
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
		BidirectedStringGraphEdge e;
		const v_idx_t v1_idx = read_1_idx;
		const v_idx_t v2_idx = read_2_idx;

		bv1.extract_seq(beg_1, end_1, e.get_seq_1_to_2());
		bv2.extract_seq(beg_2, end_2, e.get_seq_2_to_1());
		e.set_v_indices(v1_idx, v2_idx);
		e.set_dirs(dirs);

		edge_idx_t edge_idx = this->push_back_edge(e);
		_vertices[v1_idx].add_edge_idx(edge_idx);
		_vertices[v2_idx].add_edge_idx(edge_idx);
	}
};
