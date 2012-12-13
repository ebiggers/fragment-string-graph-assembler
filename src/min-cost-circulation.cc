#include "AnyStringGraph.h"

DEFINE_USAGE(
"Usage: min-cost-circulation GRAPH_FILE OUT_GRAPH_FILE\n"
"\n"
"Solves a minimum-cost circulation problem on a directed or bidirected\n"
"string graph.\n"
"\n"
"Input:\n"
"      GRAPH_FILE:   A directed or bidirected string graph in binary format.\n"
"\n"
"Output:\n"
"      OUT_GRAPH_FILE:  The resulting graph with each edge labeled with a\n"
"                       traversal count.\n"
);

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 3);
	info("Loading string graph from \"%s\"", argv[1]);
	AnyStringGraph graph(argv[1]);
	graph.min_cost_circulation();
	info("Writing string graph to \"%s\"", argv[2]);
	graph.write(argv[2]);
}
