#include "AnyStringGraph.h"

DEFINE_USAGE(
"calculate-A-statistics GRAPH_FILE OUT_GRAPH_FILE\n"
"\n"
"Calculate the A-statistic (arrival rate statistic) on each edge in\n"
"a string graph.\n"
"\n"
"Input:\n"
"      GRAPH_FILE: A string graph.\n"
"\n"
"Output:\n"
"      OUT_GRAPH_FILE:  The string graph with the A-statistics calculated.\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 3);
	const char *graph_file = argv[1];
	const char *out_graph_file = argv[2];

	info("Loading \"%s\"", graph_file);
	AnyStringGraph graph(graph_file);
	graph.calculate_A_statistics();
	info("Writing \"%s\"", out_graph_file);
	graph.write(out_graph_file);
}
