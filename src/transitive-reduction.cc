#include "Graph.h"

DEFINE_USAGE(
"Usage: transitive-reduction GRAPH_FILE OUT_GRAPH_FILE"
);

static void transitive_reduction(Graph & graph)
{
	info("Performing transitive reduction on string graph with "
	     "%zu vertices and %zu edges", graph.num_vertices(), graph.num_edges());
	info("Done performing transitive reduction");
}

int main(int argc, char **argv)
{
	USAGE_IF(argc != 3);
	info("Loading string graph from \"%s\"", argv[1]);
	Graph graph(argv[1]);
	transitive_reduction(graph);
	info("Writing string graph to \"%s\"", argv[2]);
	graph.write(argv[2]);
}
