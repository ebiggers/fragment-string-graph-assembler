#pragma once

#include "StringGraph.h"
#include "BaseVec.h"
#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <ostream>
#include <inttypes.h>

class BidirectedStringGraphEdge {
private:
	uint64_t _data;
	BaseVec _seq_1_to_2;
	BaseVec _seq_2_to_1;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _data;
		ar & _seq_1_to_2;
		ar & _seq_2_to_1;
	}
public:
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

	size_t length() const
	{
		return _seq_1_to_2.size();
	}

	unsigned long get_v1_idx() const
	{
		return (_data & 0x7fffffffULL);
	}

	unsigned long get_v2_idx() const
	{
		return (_data >> 31) & 0x7fffffffULL;
	}

	void set_v1_idx(const unsigned long v1_idx)
	{
		_data &= ~0x7fffffffULL;
		_data |= v1_idx;
	}

	void set_v2_idx(const unsigned long v2_idx)
	{
		_data &= ~(0x7fffffffULL << 31);
		_data |= (uint64_t(v2_idx) << 31);
	}

	void set_v_indices(const unsigned long v1_idx, const unsigned long v2_idx)
	{
		_data = (_data & (3ULL << 62)) | (uint64_t(v2_idx) << 31) | v1_idx;
	}

	void get_v_indices(unsigned long & v1_idx, unsigned long & v2_idx) const
	{
		v1_idx = get_v1_idx();
		v2_idx = get_v2_idx();
	}

	//DirectedStringGraphEdge(unsigned long vertex_1_id,
		  //unsigned long vertex_2_id,
		  //BaseVec seq)
		//: _vertex_1_id(vertex_1_id),
		  //_vertex_2_id(vertex_2_id),
		  //_seq(seq)
	//{ }
	//
	//
	//friend std::ostream & operator<<(std::ostream & os, const DirectedStringGraphEdge & e)
	//{
		//return os << "DirectedStringGraphEdge {_v1_idx = " << e._v1_idx << ", _v2_idx = "
			  //<< e._v2_idx << ", _seq = \"" << e._seq << "\"}";
	//}
};


class BidirectedStringGraphVertex : public StringGraphVertex {
public:
	size_t out_degree() const
	{
		unimplemented();
	}
};

class BidirectedStringGraph : public StringGraph<BidirectedStringGraphVertex,
						 BidirectedStringGraphEdge>
{
public:
	BidirectedStringGraph(size_t num_reads)
	{
		_vertices.resize(num_reads);
	}

	void write(const char *filename) const;
	void print(std::ostream & os) const;
	void print_dot(std::ostream & os) const;

	void add_edge(const unsigned long read_1_idx, const int read_1_end,
		      const unsigned long read_2_idx, const int read_2_end,
		      const BaseVec & bv,
		      const unsigned long beg, const unsigned long end)
	{
		unimplemented();
#if 0
		BidirectedStringGraphEdge e;
		unsigned long len;
		unsigned long i;
		const unsigned long v1_idx = read_1_idx;
		const unsigned long v2_idx = read_2_idx;
		e.set_v1_idx(v1_idx);
		e.set_v2_idx(v2_idx);
		BaseVec &edge_seq = e.get_seq();
		if (end > beg) {
			len = end - beg + 1;
			edge_seq.resize(len);
			for (i = 0; i < len; i++)
				edge_seq.set(i, bv[beg + i]);
		} else {
			len = beg - end + 1;
			edge_seq.resize(len);
			for (i = 0; i < len; i++)
				edge_seq.set(i, (3 ^ bv[beg - i]));
		}
		//std::cout << edge_seq << std::endl;
		unsigned long edge_idx = _edges.size();
		_edges.push_back(e);
		_vertices[v1_idx].add_edge_idx(edge_idx);
#endif
	}
};
