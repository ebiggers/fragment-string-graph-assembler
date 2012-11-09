#include "Overlap.h"
#include "util.h"

DEFINE_USAGE(
"Usage: print-overlaps OVERLAPS_FILE\n"
"\n"
"Prints the overlaps from a file.\n"
"\n"
"Input:\n"
"      OVERLAPS_FILE:  The overlaps, in binary format.\n"
);

int main(int argc, char **argv)
{
	USAGE_IF(argc != 2);
	OverlapVecVec ovv(argv[1]);
	foreach(auto overlap_set, ovv)
		foreach(const Overlap & o, overlap_set)
			std::cout << o << std::endl;
}
