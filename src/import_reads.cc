#include "BaseVecVec.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s (FASTA_FILE | FASTQ_FILE) "
					    "BASEVECVEC_FILE\n", argv[0]);
		exit(1);
	}
	info("Loading `%s'", argv[1]);
	BaseVecVec bvv(argv[1]);
	info("Loaded %zu reads", bvv.size());
	info("Writing `%s'", argv[2]);
	bvv.write(argv[2]);
}
