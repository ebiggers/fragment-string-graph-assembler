#include "BaseVecVec.h"
#include "util.h"
#include <string.h>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/serialization/vector.hpp>

#include <iostream>
#include <fstream>

const char BaseVecVec::magic[] = "BaseVecVec";
const size_t BaseVecVec::MAGIC_LEN = sizeof(magic);

const char *BaseVecVec::file_type_string(file_type ft)
{
	switch (ft) {
	case NATIVE:
		return "native BaseVecVec binary format";
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
	char buf[MAGIC_LEN];
	std::ifstream in(filename);
	in.read(buf, MAGIC_LEN);
	if (memcmp(magic, buf, MAGIC_LEN) == 0)
		return NATIVE;
	else if (buf[0] == '@')
		return FASTQ;
	else if (buf[0] == '>')
		return FASTA;
	fatal_error("`%s': Unknown file type", filename);
}

// Create a new BaseVec from the ASCII read @seq and add it to this BaseVecVec.
void BaseVecVec::push_ascii_seq(const std::string &seq)
{
	BaseVec bv;
	bv.load_from_text(seq);
	this->push_back(bv);
}

// Load the reads from a FASTA file into this BaseVecVec.
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
			seq += s;
		}
	}
	if (seq.size() != 0)
		this->push_ascii_seq(seq);
}

// Load the reads from a FASTQ file into this BaseVecVec.
void BaseVecVec::load_fastq(std::istream &in)
{
	std::string tag, seq, sep, quals;
	while (in >> tag && in >> seq && in >> sep && in >> quals) {
		boost::algorithm::trim_right(seq);
		this->push_ascii_seq(seq);
	}
}

// Load the reads from a file into this BaseVecVec.  The file may be in FASTA,
// FASTQ, or binary BaseVecVec format.
void BaseVecVec::read(const char *filename, BaseVecVec::file_type ft)
{
	if (ft == AUTODETECT)
		ft = detect_file_type(filename);
	//info("Loading \"%s\" (filetype: %s)", filename, file_type_string(ft));
	std::ifstream in(filename);
	switch (ft) {
	case NATIVE: {
			char buf[MAGIC_LEN];
			in.read(buf, MAGIC_LEN);
			boost::archive::binary_iarchive ar(in);
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
	//info("Loaded %zu reads from \"%s\")", this->size(), filename);
}

// Write the BaseVecVec to a file in FASTA, FASTQ, or native binary format.
void BaseVecVec::write(const char *filename, file_type ft) const
{
	if (ft == AUTODETECT) {
		const char *dot = strrchr(filename, '.');
		ft = NATIVE;
		if (dot) {
			if (strcasecmp(dot + 1, "fq") == 0 || strcasecmp(dot + 1, "fastq") == 0)
				ft = FASTQ;
			else if (strcasecmp(dot + 1, "fa") == 0 || strcasecmp(dot + 1, "fasta") == 0)
				ft = FASTA;
		}
	}
	//info("Writing \"%s\" [filetype: %s]", filename, file_type_string(ft));
	std::ofstream out(filename);
	switch (ft) {
	case NATIVE: {
			out.write(magic, MAGIC_LEN);
			boost::archive::binary_oarchive ar(out);
			ar << *this;
		}
		break;
	case FASTA:
		for (size_t i = 0; i < this->size(); i++) {
			const BaseVec &bv = (*this)[i];
			out << ">read_" << i + 1 << '\n';
			size_t chars_in_line = 0;
			for (size_t j = 0; j < bv.size(); j++) {
				out << BaseUtils::bin_to_ascii(bv[j]);
				if (++chars_in_line == 70) {
					out << '\n';
					chars_in_line = 0;
				}
			}
		}
		break;
	case FASTQ:
		for (size_t i = 0; i < this->size(); i++) {
			const BaseVec &bv = (*this)[i];
			out << "@read_" << i + 1 << '\n';
			for (size_t j = 0; j < bv.size(); j++)
				out << BaseUtils::bin_to_ascii(bv[j]);
			out << '\n';
			out << '+';
			out << '\n';
			for (size_t j = 0; j < bv.size(); j++)
				out << '@';
			out << '\n';
		}
		break;
	default:
		assert(0);
	}
	out.close();
	if (out.bad())
		fatal_error_with_errno("Error writing to \"%s\"", filename);
	//info("Wrote %zu reads to \"%s\"", this->size(), filename);
}
