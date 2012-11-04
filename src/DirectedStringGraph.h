#pragma once

#include "StringGraph.h"
#include "BaseVec.h"
#include <ostream>
#include <inttypes.h>

class DirectedStringGraphEdge {
private:
	unsigned _v1_idx;
	unsigned _v2_idx;
	BaseVec _seq;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _v1_idx;
		ar & _v2_idx;
		ar & _seq;
	}
public:
	BaseVec & get_seq()
	{
		return _seq;
	}

	size_t length() const
	{
		return _seq.size();
	}

	const BaseVec & get_seq() const
	{
		return _seq;
	}

	unsigned get_v1_idx() const
	{
		return _v1_idx;
	}

	unsigned get_v2_idx() const
	{
		return _v2_idx;
	}

	void set_v1_idx(const unsigned long v1_idx)
	{
		_v1_idx = v1_idx;
	}

	void set_v2_idx(const unsigned long v2_idx)
	{
		_v2_idx = v2_idx;
	}

	friend std::ostream & operator<<(std::ostream & os,
					 const DirectedStringGraphEdge & e)
	{
		return os << "DirectedStringGraphEdge {_v1_idx = " << e._v1_idx << ", _v2_idx = "
			  << e._v2_idx << ", _seq = \"" << e._seq << "\"}";
	}
};

class DirectedStringGraphVertex : public StringGraphVertex {
public:
	size_t out_degree() const
	{
		return _edge_indices.size();
	}
};

class DirectedStringGraph : public StringGraph<DirectedStringGraphVertex,
					       DirectedStringGraphEdge>
{
private:
	void add_edge(const unsigned long read_1_idx,
		      const unsigned long read_2_idx, 
		      const unsigned long dirs,
		      const BaseVec & bv,
		      const unsigned long beg, const unsigned long end)
	{
		DirectedStringGraphEdge e;
		unsigned long len;
		unsigned long i;

		//info("Add edge %lu:%c => %lu:%c (%lu, %lu)",
		     //start_read_id, (start_read_id == READ_BEGIN) ? 'B' : 'E',
		     //end_read_id, (end_read_id == READ_BEGIN) ? 'B' : 'E',
		     //beg, end);

		unsigned long v1_idx = read_1_idx * 2 + (dirs >> 1);
		unsigned long v2_idx = read_2_idx * 2 + (dirs & 1);
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
	}
public:
	DirectedStringGraph(size_t num_reads)
	{
		_vertices.resize(num_reads * 2);
	}

	DirectedStringGraph(const char *filename);

	void write(const char *filename) const;
	void print(std::ostream & os) const;
	void print_dot(std::ostream & os) const;

	void add_edge_pair(const unsigned long read_1_idx,
			   const unsigned long read_2_idx,
			   const unsigned long dirs,
			   const BaseVec & bv1,
			   const unsigned long beg_1, const unsigned long end_1,
			   const BaseVec & bv2,
			   const unsigned long beg_2, const unsigned long end_2)
	{
		add_edge(read_1_idx, read_2_idx, dirs, bv1, beg_1, end_1);
		add_edge(read_2_idx, read_1_idx, (dirs ^ 3), bv2, beg_2, end_2);
	}
};
