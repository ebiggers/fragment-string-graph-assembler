#include "BaseVecVec.h"
#include "util.h"

DEFINE_USAGE(
"Usage: convert-reads IN_READS_FILE... OUT_READS_FILE\n"
"\n"
"Converts and/or merges reads.  File formats are auto-detected.\n"
"\n"
"Input:\n"
"      IN_READS_FILE...:    One or more FASTA, FASTQ, or native binary\n"
"                           (BaseVecVec) reads files.\n"
"\n"
"Output:\n"
"      OUT_READS_FILE:      File to write the reads to.  *.fa or *.fasta\n"
"                           for FASTA, *.fq or *.fastq for FASTQ, or\n"
"                           anything else for native BaseVecVec.\n"
"\n"
"Examples:\n"
"\n"
"      Convert FASTA reads to binary reads:\n"
"            convert-reads reads.fa reads.bvv\n"
"\n"
"      Merge two FASTQ files to one binary reads file:\n"
"            convert-reads reads_1.fq reads_2.fq reads.bvv\n"
"\n"
"      Convert a binary reads file back to FASTQ format:\n"
"            convert-reads reads.bvv reads.fq\n"
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
		foreach(const BaseVec & bv, vecs[i])
			vecs[0].push_back(bv);
		vecs[i].clear();
	}

	info("Writing %zu reads to \"%s\"", vecs[0].size(), argv[num_in_files]);
	vecs[0].write(argv[num_in_files]);
}
