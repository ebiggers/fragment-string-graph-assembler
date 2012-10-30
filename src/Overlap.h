#include <boost/serialization/binary_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>
#include <vector>

class Overlap {
public:
	unsigned long read_1_idx : 24;
	unsigned long read_2_idx : 24;
	unsigned long read_1_beg : 12;
	unsigned long read_1_end : 12;
	unsigned long read_2_beg : 12;
	unsigned long read_2_end : 12;

	friend class boost::serialization::access;

	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & boost::serialization::make_binary_object(this, sizeof(*this));
	}
};

class OverlapVecVec : public std::vector<std::vector<Overlap> > {

private:
	friend class boost::serialization::access;

	template <class Archive>
	void serialize(Archive & ar, unsigned version)
	{
		ar & boost::serialization::base_object<std::vector<std::vector<Overlap> > >(*this);
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
