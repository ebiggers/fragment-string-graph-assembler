#include "BaseVecVec.h"
#include "util.h"

DEFINE_USAGE(
"Usage: convert-reads (FASTA_FILE | FASTQ_FILE | BASEVECVEC_FILE)...\n"
"                      BASEVECVEC_FILE\n"
);

int main(int argc, char *argv[])
{
	argc--;
	argv++;

	USAGE_IF(argc < 2);
	int num_in_files = argc - 1;

	std::vector<BaseVecVec> vecs(num_in_files);
	for (int i = 0; i < num_in_files; i++) {
		info("Loading reads from \"%s\"", argv[i]);
		vecs[i].read(argv[i]);
	}

	if (num_in_files > 1)
		info("Merging reads from %d files...", num_in_files);
	for (int i = 1; i < num_in_files; i++) {
		for (const BaseVec & bv : vecs[i])
			vecs[0].push_back(bv);
		vecs[i].clear();
	}

	info("Writing %zu reads to \"%s\"", vecs[0].size(), argv[num_in_files]);
	vecs[0].write(argv[num_in_files]);
}
