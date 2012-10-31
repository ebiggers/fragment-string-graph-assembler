#include <boost/serialization/binary_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <fstream>
#include <vector>
#include <set>

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
	void set(unsigned long read_1_idx, unsigned long read_1_beg,
		 unsigned long read_1_end,
		 unsigned long read_2_idx, unsigned long read_2_beg,
		 unsigned long read_2_end)
	{
		_read_1_idx = read_1_idx;
		_read_1_beg = read_1_beg;
		_read_1_end = read_1_end;
		_read_2_idx = read_2_idx;
		_read_2_beg = read_2_beg;
		_read_2_end = read_2_end;
	}

	friend std::ostream & operator<<(std::ostream & os, const Overlap & o)
	{
		os << "Overlap : Read " << o._read_1_idx << " [" << o._read_1_beg
		   << ", " << o._read_1_end << "], Read " << o._read_2_idx
		   << " [" << o._read_1_beg << ", " << o._read_2_end << "]";
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
private:
	typedef std::vector<std::set<Overlap> > BaseT;

	friend class boost::serialization::access;

	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & boost::serialization::base_object<BaseT>(*this);
	}
public:
	OverlapVecVec() { }

	OverlapVecVec(const char *filename) {
		info("Loading overlaps from \"%s\"", filename);
		std::ifstream in(filename);
		boost::archive::binary_iarchive ar(in);
		ar >> *this;
		info("Done loading overlaps from \"%s\"", filename);
	}

	void write(const char *filename) {
		info("Writing overlaps to \"%s\"", filename);
		std::ofstream out(filename);
		boost::archive::binary_oarchive ar(out);
		ar << *this;
		out.close();
		if (out.bad())
			fatal_error_with_errno("Error writing to \"%s\"", filename);
		info("Done writing overlaps to \"%s\"", filename);
	}
};
