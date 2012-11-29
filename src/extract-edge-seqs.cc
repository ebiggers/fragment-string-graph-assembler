#include "AnyStringGraph.h"
#include "util.h"

DEFINE_USAGE(
"Usage: extract-edge-seqs GRAPH_FILE OUT_CONTIGS_FILE\n"
"\n"
"Extracts the sequences from the edges of a string graph.\n"
"\n"
"Input:\n"
"      GRAPH_FILE:  The string graph.\n"
"      OUT_CONTIGS_FILE:  A file to which the edge sequences are to be\n"
"                         extracted.\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 3);
	const char *graph_file = argv[1];
	const char *out_contigs_file = argv[2];

	AnyStringGraph graph(graph_file);

	BaseVecVec bvv;

	graph.extract_edge_seqs(bvv);
	info("Writing %zu contigs to \"%s\"", bvv.size(), out_contigs_file);
	bvv.write(out_contigs_file);
}
