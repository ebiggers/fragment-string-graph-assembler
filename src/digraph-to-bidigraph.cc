#include "BidirectedStringGraph.h"
#include "DirectedStringGraph.h"
#include "util.h"

DEFINE_USAGE(
"Usage: digraph-to-bidigraph DIGRAPH_FILE OUT_BIDIGRAPH_FILE\n"
"\n"
"Turns a directed string graph into a bidirected string graph.\n"
"\n"
"Input:\n"
"      DIGRAPH_FILE:   A directed string graph in binary format.\n"
"\n"
"Output:\n"
"      OUT_BIDIGRAPH_FILE:  A bidirected string graph in binary format.\n"
);

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 3);
	info("Loading directed string graph from \"%s\"", argv[1]);
	DirectedStringGraph digraph(argv[1]);

	BidirectedStringGraph bidigraph(digraph.num_vertices() / 2);

	bidigraph.build_from_digraph(digraph);

	info("Writing bidirected string graph to \"%s\"", argv[2]);
	bidigraph.write(argv[2]);
}
