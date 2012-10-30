#include "BaseVecVec.h"

const char * const BaseVecVec::magic = "BaseVecVec";

const char *BaseVecVec::file_type_string(file_type ft)
{
	switch (ft) {
	case NATIVE:
		return "native (BaseVecVec binary format)";
	case FASTA:
		return "FASTA";
	case FASTQ:
		return "FASTQ";
	default:
		return "Unknown";
	}
}

BaseVecVec::file_type BaseVecVec::detect_file_type(const char *filename)
{
	char buf[sizeof(magic)];
	std::ifstream in(filename);
	in.read(buf, sizeof(buf));
	if (memcmp(magic, buf, sizeof(buf)) == 0) {
		return NATIVE;
	} else if (buf[0] == '@') {
		return FASTQ;
	} else if (buf[0] == '>') {
		return FASTA;
	}
	fatal_error("`%s': Unknown file type", magic);
}

void BaseVecVec::push_ascii_seq(const std::string &seq)
{
	info("Push back seq %s", seq.c_str());
	BaseVec bv;
	bv.load_from_text(seq);
	this->push_back(bv);
	bv.leak();
}

void BaseVecVec::load_fasta(std::istream &in)
{
	std::string seq;
	std::string s;
	BaseVec cur_bv;
	while (in >> s) {
		boost::algorithm::trim_right(s);
		if (s[0] == '>') {
			if (seq.size() != 0) {
				this->push_ascii_seq(seq);
				seq.clear();
			}
		} else {
			seq.insert(seq.end(), s.begin(), s.end());
		}
	}
	if (seq.size() != 0) {
		this->push_ascii_seq(seq);
	}
}

void BaseVecVec::load_fastq(std::istream &in)
{
	std::string tag, seq, sep, quals;
	while (in >> tag && in >> seq && in >> sep && in >> quals) {
		boost::algorithm::trim_right(seq);
		this->push_ascii_seq(seq);
	}
}

BaseVecVec::BaseVecVec(const char *filename, file_type ft)
{
	if (ft == AUTODETECT)
		ft = detect_file_type(filename);
	info("Loading `%s' [filetype: %s]", filename, file_type_string(ft));
	std::ifstream in(filename);
	switch (ft) {
	case NATIVE: {
			boost::archive::text_iarchive ar(in);
			ar >> *this;
		}
		break;
	case FASTA:
		load_fasta(in);
		break;
	case FASTQ:
		load_fastq(in);
		break;
	default:
		assert(0);
	}
	info("Loaded `%s' (%zu reads)", filename, size());
}

void BaseVecVec::write(const char *filename, file_type ft)
{
	info("Writing `%s' [filetype: %s]", filename, file_type_string(ft));
	switch (ft) {
	case NATIVE: {
			std::ofstream out(filename);
			boost::archive::text_oarchive ar(out);
			ar << *this;
		}
	case FASTA:
	case FASTQ:
		unimplemented();
	default:
		assert(0);
	}
}
