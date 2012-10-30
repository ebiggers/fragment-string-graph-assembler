#include "BaseVecVec.h"
#include "Overlap.h"

static void compute_overlaps(const BaseVecVec &bvv, OverlapVecVec &ovv)
{
	ovv.clear();
	ovv.resize(bvv.size());
	for (size_t i = 0; i < bvv.size(); i++) {
		Overlap o;
		o.read_1_idx = 0;
		o.read_2_idx = 1;
		o.read_1_beg = 2;
		o.read_2_beg = 3;
		o.read_1_end = 4;
		o.read_2_end = 5;
		ovv[i].push_back(o);
	}
}

int main(int argc, char *argv[])
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s BVV_FILE OVERLAPS_FILE\n", argv[0]);
		exit(1);
	}
	BaseVecVec bvv(argv[1]);
	OverlapVecVec ovv;
	compute_overlaps(bvv, ovv);
	ovv.write(argv[2]);
}
