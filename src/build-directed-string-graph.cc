#include "Overlap.h"
#include "BaseVecVec.h"
#include "DirectedStringGraph.h"

DEFINE_USAGE(
"Usage: build-directed-string-graph READS_FILE OVERLAPS_FILE DIGRAPH_FILE\n"
"\n"
"Builds a directed string graph.\n"
"\n"
"Input:\n"
"      READS_FILE:      The set of reads from which the overlaps were\n"
"                       computed.\n"
"      OVERLAPS_FILE:   The overlaps between the reads.\n"
"\n"
"Output:\n"
"      DIGRAPH_FILE:    File containing the directed string graph\n"
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

	DirectedStringGraph graph(bvv.size());

	info("Building directed string graph from overlaps");
	graph.build(bvv, ovv);

	info("Writing directed string graph to \"%s\"", graph_file);
	graph.write(graph_file);

	assert(ovv.size() == bvv.size());

	info("Done");
}
