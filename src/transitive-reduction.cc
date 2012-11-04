#include "DirectedStringGraph.h"
#include <algorithm>

DEFINE_USAGE(
"Usage: transitive-reduction GRAPH_FILE OUT_GRAPH_FILE\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 3);
	info("Loading directed string graph from \"%s\"", argv[1]);
	DirectedStringGraph graph(argv[1]);
	graph.transitive_reduction();
	info("Writing directed string graph to \"%s\"", argv[2]);
	graph.write(argv[2]);
}
