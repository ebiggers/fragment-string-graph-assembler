#include "BaseVec.h"
#include "util.h"

#include <string.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_iarchive.hpp>

class BaseVecVec : std::vector<BaseVec> {
public:
	enum file_type {
		NATIVE,
		FASTA,
		FASTQ,
		AUTODETECT,
	};
private:
	static const char * const magic;

	file_type detect_file_type(const char *filename) {
		char buf[sizeof(magic)];
		std::ifstream in(filename);
		in.read(buf, sizeof(buf));
		if (memcmp(magic, buf, sizeof(magic)) == 0) {
			return NATIVE;
		} else if (buf[0] == '@') {
			return FASTQ;
		} else if (buf[0] == '>') {
			return FASTA;
		}
		fatal_error("`%s': Unknown file type", magic);
	}
public:
	BaseVecVec(const char *filename, file_type ft = AUTODETECT) {
		if (ft == AUTODETECT)
			ft = detect_file_type(filename);
		std::ifstream in(filename);
		boost::archive::text_iarchive ar(in);
	}
};
