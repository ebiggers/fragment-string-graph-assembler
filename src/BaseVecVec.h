#include "BaseVec.h"
#include "util.h"

#include <string.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <string>

class BaseVecVec : public std::vector<BaseVec> {
public:
	enum file_type {
		NATIVE,
		FASTA,
		FASTQ,
		AUTODETECT,
	};
private:
	static const char magic[];
	static const size_t MAGIC_LEN;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version) {
		ar & boost::serialization::base_object<std::vector<BaseVec> >(*this);
	}

	static const char *file_type_string(file_type ft);
	static file_type detect_file_type(const char *filename);
	void push_ascii_seq(const std::string &seq);
	void load_fasta(std::istream &in);
	void load_fastq(std::istream &in);

public:
	BaseVecVec(const char *filename, file_type ft = AUTODETECT);
	void write(const char *filename, file_type ft = NATIVE);
};
