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

class Overlap {
private:
	unsigned long _read_1_idx : 24;
	unsigned long _read_1_beg : 12;
	unsigned long _read_1_end : 12;
	unsigned long _read_2_idx : 24;
	unsigned long _read_2_beg : 12;
	unsigned long _read_2_end : 12;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & boost::serialization::make_binary_object(this, sizeof(*this));
	}

public:
	void set(unsigned long read_1_idx,
		 unsigned long read_1_beg,
		 unsigned long read_1_end,
		 unsigned long read_2_idx,
		 unsigned long read_2_beg,
		 unsigned long read_2_end)
	{
		assert(read_1_end > read_1_beg ||
		       read_2_end > read_2_beg);

		_read_1_idx = read_1_idx;
		_read_1_beg = read_1_beg;
		_read_1_end = read_1_end;
		_read_2_idx = read_2_idx;
		_read_2_beg = read_2_beg;
		_read_2_end = read_2_end;
	}

	void set_indices(unsigned long read_1_idx,
			 unsigned long read_2_idx)
	{
		_read_1_idx = read_1_idx;
		_read_2_idx = read_2_idx;
	}

	void get(unsigned long & read_1_idx,
		 unsigned long & read_1_beg,
		 unsigned long & read_1_end,
		 unsigned long & read_2_idx,
		 unsigned long & read_2_beg,
		 unsigned long & read_2_end) const
	{
		read_1_idx = _read_1_idx;
		read_1_beg = _read_1_beg;
		read_1_end = _read_1_end;
		read_2_idx = _read_2_idx;
		read_2_beg = _read_2_beg;
		read_2_end = _read_2_end;
	}

	void get_indices(unsigned long & read_1_idx,
			 unsigned long & read_2_idx) const
	{
		read_1_idx = _read_1_idx;
		read_2_idx = _read_2_idx;
	}

	friend std::ostream & operator<<(std::ostream & os, const Overlap & o)
	{
		os << "Overlap { Read " << o._read_1_idx << ": [" << o._read_1_beg
		   << ", " << o._read_1_end << "], Read " << o._read_2_idx
		   << ": [" << o._read_2_beg << ", " << o._read_2_end << "] }";
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

	OverlapVecVec(const char *filename)
	{
		std::ifstream in(filename);
		boost::archive::binary_iarchive ar(in);
		ar >> *this;
	}

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
			      const unsigned pos1,
			      const unsigned pos2,
			      const unsigned len,
			      const bool is_rc1, const bool is_rc2,
			      const char *description = "SEED");

extern void assert_overlap_valid(const Overlap & o, const BaseVecVec & bvv,
				 const unsigned min_overlap_len,
				 const unsigned max_edits);
