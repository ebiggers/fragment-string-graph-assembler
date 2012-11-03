#include "Graph.h"
#include <fstream>
#include <ostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

Graph::Graph(const char *filename)
{
	std::ifstream in(filename);
	boost::archive::binary_iarchive ar(in);
	ar >> *this;
}

void Graph::write(const char *filename) const
{
	std::ofstream out(filename);
	boost::archive::binary_oarchive ar(out);
	ar << *this;
	out.close();
	if (out.bad())
		fatal_error_with_errno("Error writing to \"%s\"", filename);
}

void Graph::print(std::ostream & os) const
{
	for (size_t i = 0; i < _vertices.size(); i++) {
		//os << "[Vertex " << i << "]" << '\n';
		const GraphVertex & v = _vertices[i];
		const std::vector<unsigned long> & edge_indices = v.edge_indices();
		for (size_t j = 0; j < edge_indices.size(); j++) {
			const GraphEdge & e = _edges[edge_indices[j]];
			unsigned v1_idx = e.get_v1_idx();
			size_t read_1_idx = v1_idx / 2 + 1;
			char read_1_dir = (v1_idx & 1) ? 'E' : 'B';
			unsigned v2_idx = e.get_v2_idx();
			size_t read_2_idx = v2_idx / 2 + 1;
			char read_2_dir = (v2_idx & 1) ? 'E' : 'B';
			os << read_1_idx << '.' << read_1_dir << " -> "
			   << read_2_idx << '.' << read_2_dir
			   << '\t' << e.get_seq() << '\n';
		}
	}
	os << std::flush;
}

void Graph::print_dot(std::ostream & os) const
{
	os << "digraph {\n";
	os << "\tnode [shape = oval];\n";
	for (size_t i = 0; i < _vertices.size(); i++) {
		size_t read_idx = i / 2;
		char read_dir = (i & 1) ? 'E' : 'B';
		os << "\tv" << i << " [label = \"" << (read_idx + 1) << '.' << read_dir << "\"];\n";
	}
	for (size_t i = 0; i < _edges.size(); i++) {
		const GraphEdge & e = _edges[i];
		os << "\tv" << e.get_v1_idx() << " -> " << "v" << e.get_v2_idx()
		   << " [ label = \"" << e.get_seq().size() << "\" ];\n";
	}
	os << "}" << std::endl;
}
