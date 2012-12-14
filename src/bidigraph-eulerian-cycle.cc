#include "BidirectedStringGraph.h"

DEFINE_USAGE(
"Usage: bidigraph-eulerian-cycle BIDIGRAPH_FILE OUT_CYCLE_FILE\n"
"\n"
"Finds an Eulerian cycle in a bidirected graph.\n"
"\n"
"Input:\n"
"      BIDIGRAPH_FILE:   A bidirected string graph in binary format.\n"
"\n"
"Output:\n"
"      OUT_CYCLE_FILE:  The resulting Eulerian cycle as a vector of edge\n"
"                       indices (in binary format)\n"
);

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 3);
	info("Loading bidirected graph graph from \"%s\"", argv[1]);
	BidirectedStringGraph graph(argv[1]);

	std::vector<size_t> cycle;

	graph.eulerian_cycle(cycle);

	info("Writing Eulerian cycle to \"%s\"", argv[2]);
	{
		std::ofstream out(argv[2]);
		boost::archive::binary_oarchive ar(out);
		ar << cycle;
		out.close();
		if (!out)
			fatal_error_with_errno("Error writing to \"%s\"", argv[2]);
	}
}
