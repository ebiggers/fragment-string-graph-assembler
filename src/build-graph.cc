#include "Overlap.h"
#include "util.h"

DEFINE_USAGE(
"Usage: build-graph OVERLAPS_FILE GRAPH_FILE\n"
);

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 3);
	const char *overlaps_file = argv[1];
	const char *graph_file = argv[2];

	info("Loading overlaps from \"%s\"", overlaps_file);
	OverlapVecVec ovv(overlaps_file);
}
