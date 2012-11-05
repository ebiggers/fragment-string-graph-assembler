#include "AnyStringGraph.h"

DEFINE_USAGE(
"Usage: collapse-unbranched-paths GRAPH_FILE OUT_GRAPH_FILE\n"
"\n"
"Collapse unbranched paths in a string graph.\n"
"\n"
"Input:\n"
"      GRAPH_FILE:   A directed or bidirected string graph in binary format.\n"
"\n"
"Output:\n"
"      OUT_GRAPH_FILE:  The graph with unbranched paths collapsed.\n"
);

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 3);
	info("Loading string graph from \"%s\"", argv[1]);
	AnyStringGraph graph(argv[1]);
	graph.collapse_unbranched_paths();
	info("Writing string graph to \"%s\"", argv[2]);
	graph.write(argv[2]);
}
