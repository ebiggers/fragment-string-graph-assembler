#include "Graph.h"

DEFINE_USAGE(
"Usage: print-graph BIN_GRAPH_FILE\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 2);
	Graph(argv[1]).print();
}
