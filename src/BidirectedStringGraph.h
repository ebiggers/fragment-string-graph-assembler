#pragma once

#include "StringGraph.h"
#include "BaseVec.h"
#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <ostream>
#include <inttypes.h>

class BidirectedStringGraphEdge : public StringGraphEdge {
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

	size_t length() const
	{
		return _seq_1_to_2.size();
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

	void set_dirs(unsigned dirs)
	{
		_data = (_data & ~(3ULL << 62)) | (uint64_t(dirs) << 62);
	}

	//friend std::ostream & operator<<(std::ostream & os,
					 //const BidirectedStringGraphEdge & e)
	//{
		//unimplemented();
	//}
	void print_dot(std::ostream & os, size_t e_idx) const
	{
		unimplemented();
		//os << "v" << get_v1_idx() << " -> "
		   //<< "v" << get_v2_idx()
		   //<< " [ label = \"" << length() << "\" ]";
	}
};


class BidirectedStringGraphVertex : public StringGraphVertex {
public:
	size_t out_degree() const 
	{
		unimplemented();
	}

	void print_dot(std::ostream & os, size_t v_idx) const
	{
		os << "v" << v_idx << " [label = \"" << (v_idx + 1) << "\"]";
	}
};

class BidirectedStringGraph : public StringGraph<BidirectedStringGraphVertex,
						 BidirectedStringGraphEdge,
						 BidirectedStringGraph>
{
public:
	BidirectedStringGraph(size_t num_reads)
	{
		if (!enough_v_indices(num_reads * 2))
			fatal_error("Too many reads (%zu)", num_reads);
		_vertices.resize(num_reads);
	}

	BidirectedStringGraph(const char *filename)
	{
		this->read(filename);
	}

	void transitive_reduction();

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
