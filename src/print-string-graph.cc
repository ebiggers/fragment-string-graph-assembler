#include "AnyStringGraph.h"
#include <iostream>
#include <getopt.h>

DEFINE_USAGE(
"Usage: print-string-graph [--dot] GRAPH_FILE\n"
"\n"
"Prints a directed or bidirected string graph.\n"
"\n"
"Input:\n"
"      GRAPH_FILE:  A directed or bidirected string graph in binary format.\n"
"\n"
"Options:\n"
"   --dot    Print the graph in DOT format.\n"
"   --stats  Print statistics about the graph.\n"
);

static const char *optstring = "dsh";
static const struct option longopts[] = {
	{"dot", no_argument, NULL, 'd'},
	{"stats", no_argument, NULL, 's'},
	END_LONGOPTS
};

int main(int argc, char **argv)
{
	int c;
	bool dot = false;
	bool stats = false;
	for_opt(c) {
		switch (c) {
		case 'd':
			dot = true;
			break;
		case 's':
			stats = true;
			break;
		PROCESS_OTHER_OPTS
		}
	}
	argc -= optind;
	argv += optind;
	USAGE_IF(argc != 1);
	AnyStringGraph graph(argv[0]);
	if (dot)
		graph.print_dot(std::cout);
	else if (stats)
		graph.print_stats(std::cout);
	else
		graph.print(std::cout);
}
