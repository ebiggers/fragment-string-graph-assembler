#include "BidirectedStringGraph.h"
#include "DirectedStringGraph.h"
#include "util.h"

DEFINE_USAGE(
"Usage: bidigraph-to-digraph BIDIGRAPH_FILE OUT_DIGRAPH_FILE\n"
"\n"
"Turns a bidirected string graph into a directed string graph.\n"
"\n"
"Input:\n"
"      BIDIGRAPH_FILE:   A bidirected string graph in binary format.\n"
"\n"
"Output:\n"
"      OUT_DIGRAPH_FILE:  A directed string graph in binary format.\n"
);

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 3);
	info("Loading bidirected string graph from \"%s\"", argv[1]);
	BidirectedStringGraph bidigraph(argv[1]);

	DirectedStringGraph digraph(bidigraph.num_vertices());

	digraph.build_from_bidigraph(bidigraph);

	info("Writing directed string graph to \"%s\"", argv[2]);
	digraph.write(argv[2]);
}
