#include "Overlap.h"
#include "util.h"

DEFINE_USAGE(
"Usage: print-overlaps BIN_OVERLAPS_FILE\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 2);
	OverlapVecVec ovv(argv[1]);
	for (auto overlap_set : ovv)
		for (const Overlap & o : overlap_set)
			std::cout << o << std::endl;
}
