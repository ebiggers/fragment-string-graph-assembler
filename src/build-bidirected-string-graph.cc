#include "Overlap.h"
#include "BaseVecVec.h"
#include "BidirectedStringGraph.h"

DEFINE_USAGE(
"Usage: build-bidirected-string-graph READS_FILE OVERLAPS_FILE BIDIGRAPH_FILE\n"
"\n"
"Builds a bidirected string graph.\n"
"\n"
"Input:\n"
"      READS_FILE:      The set of reads from which the overlaps were\n"
"                       computed.\n"
"      OVERLAPS_FILE:   The overlaps between the reads.\n"
"\n"
"Output:\n"
"      BIDIGRAPH_FILE:  File containing the bidirected string graph\n"
"                       in binary format.\n"
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

	BidirectedStringGraph graph(bvv.size());

	info("Building bidirected string graph from overlaps");
	graph.build(bvv, ovv);

	info("Writing bidirected string graph to \"%s\"", graph_file);
	graph.write(graph_file);

	assert(ovv.size() == bvv.size());

	info("Done");
}
