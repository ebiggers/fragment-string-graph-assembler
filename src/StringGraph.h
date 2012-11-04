#pragma once

#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "util.h"

class BaseVec;

class StringGraphVertex {
protected:
	std::vector<unsigned long> _edge_indices;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _edge_indices;
	}

	StringGraphVertex() { }
public:
	void add_edge_idx(const unsigned long edge_idx)
	{
		_edge_indices.push_back(edge_idx);
	}

	const std::vector<unsigned long> & edge_indices() const
	{
		return _edge_indices;
	}

	std::vector<unsigned long> & edge_indices()
	{
		return _edge_indices;
	}

	friend std::ostream & operator<<(std::ostream & os,
					 const StringGraphVertex & v)
	{
		os << "StringGraphVertex {_edge_indices = [";
		for (unsigned long idx : v.edge_indices())
			os << idx << ", ";
		return os << "] }";
	}
};

template<class VERTEX_t, class EDGE_t>
class StringGraph {
protected:
	std::vector<VERTEX_t> _vertices;
	std::vector<EDGE_t> _edges;

	friend class boost::serialization::access;

	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & _vertices;
		ar & _edges;
	}

	StringGraph() { }
public:

	static const unsigned long TAG_F_B = 0x0;
	static const unsigned long TAG_F_E = 0x2;
	static const unsigned long TAG_G_B = 0x0;
	static const unsigned long TAG_G_E = 0x1;

	StringGraph(const char *filename)
	{
		unimplemented();
	}

	enum {
		READ_BEGIN = 0,
		READ_END = 1
	};

	std::vector<EDGE_t> & edges()
	{
		return _edges;
	}

	std::vector<VERTEX_t> & vertices()
	{
		return _vertices;
	}

	size_t num_edges() const
	{
		return _edges.size();
	}

	size_t num_vertices() const
	{
		return _vertices.size();
	}

	void write(const char *filename) const
	{
		unimplemented();
	}

	void print(std::ostream & os) const
	{
		for (const VERTEX_t & v : _vertices) {
			for (const unsigned long edge_idx : v.edge_indices()) {
				os << _edges[edge_idx] << '\n';
			}
		}
		os << std::flush;
	}

	void print_dot(std::ostream & os) const
	{
		os << "digraph {\n";
		os << "\tnode [shape = oval];\n";
		for (size_t i = 0; i < _vertices.size(); i++) {
			os << "\t";
			_vertices[i].print_dot(os, i);
			os << ";\n";
		}

		for (size_t i = 0; i < _edges.size(); i++) {
			os << "\t";
			_edges[i].print_dot(os, i);
			os << ";\n";
		}
		os << "}" << std::endl;
	}

	void add_edge_pair(const unsigned long read_1_idx,
			   const unsigned long read_2_idx,
			   const unsigned long dirs,
			   const BaseVec & bv1,
			   const unsigned long beg_1, const unsigned long end_1,
			   const BaseVec & bv2,
			   const unsigned long beg_2, const unsigned long end_2)
	{
		unimplemented();
	}

	void transitive_reduction()
	{
		unimplemented();
	}
};
