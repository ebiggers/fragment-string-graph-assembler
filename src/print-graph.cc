#include "Graph.h"
#include <iostream>
#include <getopt.h>

DEFINE_USAGE(
"Usage: print-graph [--dot] BIN_GRAPH_FILE\n"
);

static const char *optstring = "dh";
static const struct option longopts[] = {
	{"dot", no_argument, NULL, 'd'},
	END_LONGOPTS
};

int main(int argc, char **argv)
{
	int c;
	bool dot = false;
	for_opt(c) {
		switch (c) {
		case 'd':
			dot = true;
			break;
		PROCESS_OTHER_OPTS
		}
	}
	argc -= optind;
	argv += optind;
	USAGE_IF(argc != 1);
	DirectedStringGraph graph(argv[0]);
	if (dot)
		graph.print_dot(std::cout);
	else
		graph.print(std::cout);
}
