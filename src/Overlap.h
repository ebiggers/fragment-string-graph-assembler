#pragma once

#include <boost/serialization/binary_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <fstream>
#include <vector>
#include <set>
#include "util.h"
#include <assert.h>

class BaseVecVec;
class BaseVec;

//
// Represents an overlap between two reads:
//
// The bases in the read at _read_1_idx, beginning at _read_1_beg and ending at
// _read_2_end (both inclusive), match the bases in the read at _read_2_idx,
// beginning at _read_2_beg and ending at _read_2_end (both inclusive).  If _rc
// is 1, it is actually the reverse-complement sequence that is matched.
//
class Overlap {
private:
	unsigned long _read_1_idx : 24;
	unsigned long _read_1_beg : 12;
	unsigned long _read_1_end : 12;
	unsigned long _read_2_idx : 24;
	unsigned long _read_2_beg : 12;
	unsigned long _read_2_end : 12;
	unsigned long _rc         : 1;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & boost::serialization::make_binary_object(this, sizeof(*this));
	}

public:
	static const size_t MAX_READ_IDX = (1 << 24) - 1;
	static const size_t MAX_READ_LEN = (1 << 12) - 1;

	typedef unsigned int read_idx_t;
	typedef unsigned int read_pos_t;

	void set(const read_idx_t read_1_idx,
		 const read_pos_t read_1_beg,
		 const read_pos_t read_1_end,
		 const read_idx_t read_2_idx,
		 const read_pos_t read_2_beg,
		 const read_pos_t read_2_end,
		 const bool rc)
	{
		assert(read_1_end >= read_1_beg);
		assert(read_2_end >= read_2_end);
		assert(read_1_beg < MAX_READ_LEN);
		assert(read_1_end < MAX_READ_LEN);
		assert(read_2_beg < MAX_READ_LEN);
		assert(read_2_end < MAX_READ_LEN);
		assert(read_1_idx < MAX_READ_IDX);
		assert(read_2_idx < MAX_READ_IDX);

		_read_1_idx = read_1_idx;
		_read_1_beg = read_1_beg;
		_read_1_end = read_1_end;
		_read_2_idx = read_2_idx;
		_read_2_beg = read_2_beg;
		_read_2_end = read_2_end;
		_rc         = rc;
	}

	void set_indices(read_idx_t read_1_idx,
			 read_idx_t read_2_idx)
	{
		_read_1_idx = read_1_idx;
		_read_2_idx = read_2_idx;
	}

	void get(read_idx_t & read_1_idx,
		 read_pos_t & read_1_beg,
		 read_pos_t & read_1_end,
		 read_idx_t & read_2_idx,
		 read_pos_t & read_2_beg,
		 read_pos_t & read_2_end,
		 bool       & rc) const
	{
		read_1_idx = _read_1_idx;
		read_1_beg = _read_1_beg;
		read_1_end = _read_1_end;
		read_2_idx = _read_2_idx;
		read_2_beg = _read_2_beg;
		read_2_end = _read_2_end;
		rc = _rc;
	}

	void get_indices(read_idx_t & read_1_idx,
			 read_idx_t & read_2_idx) const
	{
		read_1_idx = _read_1_idx;
		read_2_idx = _read_2_idx;
	}

	friend std::ostream & operator<<(std::ostream & os, const Overlap & o)
	{
		os << "Overlap { Read " << (o._read_1_idx + 1) << ": [" << o._read_1_beg
		   << ", " << o._read_1_end << "], Read " << (o._read_2_idx + 1)
		   << ": [" << o._read_2_beg << ", " << o._read_2_end
		   << "], rc = " << o._rc << " }";
		return os;
	}

	friend bool operator==(const Overlap & o1, const Overlap & o2)
	{
		return memcmp(&o1, &o2, sizeof(Overlap)) == 0;
	}

	friend bool operator<(const Overlap & o1, const Overlap & o2)
	{
		return memcmp(&o1, &o2, sizeof(Overlap)) < 0;
	}
};

//
// A set of overlaps for each read.
//
class OverlapVecVec : public std::vector<std::set<Overlap> > {
public:
	typedef std::set<Overlap> OverlapSet;
private:
	typedef std::vector<OverlapSet > BaseT;

	friend class boost::serialization::access;

	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & boost::serialization::base_object<BaseT>(*this);
	}
public:
	OverlapVecVec() { }

	// Read the overlaps from a file.
	OverlapVecVec(const char *filename)
	{
		std::ifstream in(filename);
		boost::archive::binary_iarchive ar(in);
		ar >> *this;
	}

	// Write the overlaps to a file.
	void write(const char *filename)
	{
		std::ofstream out(filename);
		boost::archive::binary_oarchive ar(out);
		ar << *this;
		out.close();
		if (out.bad())
			fatal_error_with_errno("Error writing to \"%s\"", filename);
	}
};

extern void assert_seed_valid(const BaseVec & bv1,
			      const BaseVec & bv2,
			      const Overlap::read_pos_t pos1,
			      const Overlap::read_pos_t pos2,
			      const Overlap::read_pos_t len,
			      const bool is_rc1, const bool is_rc2,
			      const char *description = "SEED");

extern void assert_overlap_valid(const Overlap & o, const BaseVecVec & bvv,
				 const unsigned min_overlap_len,
				 const unsigned max_edits);
