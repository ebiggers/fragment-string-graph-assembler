#include <vector>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include "util.h"

class GraphEdge {
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

	const BaseVec & get_seq() const
	{
		return _seq;
	}

	void set_v1_idx(const unsigned long v1_idx)
	{
		_v1_idx = v1_idx;
	}

	void set_v2_idx(const unsigned long v2_idx)
	{
		_v2_idx = v2_idx;
	}
	//GraphEdge(unsigned long vertex_1_id,
		  //unsigned long vertex_2_id,
		  //BaseVec seq)
		//: _vertex_1_id(vertex_1_id),
		  //_vertex_2_id(vertex_2_id),
		  //_seq(seq)
	//{ }
};

class GraphVertex {
private:
	std::vector<unsigned long> _edge_indices;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _edge_indices;
	}
public:
	void add_edge_idx(const unsigned long edge_idx)
	{
		_edge_indices.push_back(edge_idx);
	}

	size_t out_degree() const {
		return _edge_indices.size();
	}
};

class Graph {
private:
	std::vector<GraphVertex> _vertices;
	std::vector<GraphEdge> _edges;

	friend class boost::serialization::access;

	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _vertices;
		ar & _edges;
	}

public:
	Graph(size_t num_reads)
	{
		_vertices.resize(num_reads * 2);
	}

	Graph(const char *filename)
	{
		std::ifstream in(filename);
		boost::archive::binary_iarchive ar(in);
		ar >> *this;
	}

	enum {
		READ_BEGIN = 0,
		READ_END = 1
	};


	void write(const char *filename) const
	{
		std::ofstream out(filename);
		boost::archive::binary_oarchive ar(out);
		ar << *this;
		out.close();
		if (out.bad())
			fatal_error_with_errno("Error writing to \"%s\"", filename);
	}

	void add_edge(const unsigned long start_read_id, const int start_read_end,
		      const unsigned long end_read_id, const int end_read_end,
		      const BaseVec & bv,
		      const unsigned long beg, const unsigned long end)
	{
		GraphEdge e;
		unsigned long len;
		unsigned long i;

		//info("Add edge %lu:%c => %lu:%c (%lu, %lu)",
		     //start_read_id, (start_read_id == READ_BEGIN) ? 'B' : 'E',
		     //end_read_id, (end_read_id == READ_BEGIN) ? 'B' : 'E',
		     //beg, end);

		unsigned long v1_idx = start_read_id * 2 + start_read_end;
		unsigned long v2_idx = end_read_id * 2 + end_read_end;
		e.set_v1_idx(v1_idx);
		e.set_v2_idx(v2_idx);
		BaseVec &edge_seq = e.get_seq();
		if (end > beg) {
			len = end - beg + 1;
			edge_seq.resize(len);
			for (i = 0; i < len; i++)
				edge_seq.set(i, bv[i + beg]);
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
};
