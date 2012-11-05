#include "Overlap.h"
#include "BaseVecVec.h"
#include "DirectedStringGraph.h"

DEFINE_USAGE(
"Usage: build-directed-string-graph READS_FILE OVERLAPS_FILE GRAPH_FILE\n"
);

int main(int argc, char *argv[])
{
	USAGE_IF(argc != 4);
	const char *reads_file = argv[1];
	const char *overlaps_file = argv[2];
	const char *graph_file = argv[3];

	info("Reading reads from \"%s\"", reads_file);
	BaseVecVec bvv(reads_file);
	info("Loaded %zu reads from \"%s\"", bvv.size(), reads_file);

	info("Loading overlaps from \"%s\"", overlaps_file);
	OverlapVecVec ovv(overlaps_file);
	info("Done loading overlaps");

	DirectedStringGraph graph(bvv.size());

	info("Building directed string graph from overlaps");
	graph.build(bvv, ovv);

	info("Writing directed string graph to \"%s\"", graph_file);
	graph.write(graph_file);

	assert(ovv.size() == bvv.size());

	info("Done");
}
