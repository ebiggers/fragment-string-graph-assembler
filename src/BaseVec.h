#pragma once

#include "BaseUtils.h"
#include "util.h"

#include <boost/serialization/split_member.hpp>
#include <string>
#include <ostream>
#include <assert.h>

//
// Vector of DNA bases, stored in binary format (2 bits per base).
//
class BaseVec {
public:
	typedef unsigned char storage_type;
	typedef unsigned size_type;
private:
	size_type _size;
	storage_type *_bases;
	static const size_type BITS_PER_BASE = 2;
	static const size_type BASES_PER_BYTE = 8 / BITS_PER_BASE;
	static const size_type BASES_PER_STORAGE_TYPE =
				sizeof(storage_type) * BASES_PER_BYTE;
	static const storage_type BASE_MASK =
			(static_cast<storage_type>(1) << BITS_PER_BASE) - 1;

public:

	// Return the number of bases in this BaseVec.
	size_type size() const
	{
		return _size;
	}

	// Get the binary base in this BaseVec at index @idx.
	unsigned char operator[](size_type idx) const
	{
		assert2(idx < _size);
		size_type slot = idx / BASES_PER_STORAGE_TYPE;
		size_type offset = (idx % BASES_PER_STORAGE_TYPE) * BITS_PER_BASE;
		return (_bases[slot] >> offset) & BASE_MASK;
	}

	// Set base @idx in this BaseVec to the binary base @base.
	void set(size_type idx, unsigned char base)
	{
		assert2(base < 4);
		assert2(idx < _size);
		size_type slot = idx / BASES_PER_STORAGE_TYPE;
		size_type offset = (idx % BASES_PER_STORAGE_TYPE) *
				    BITS_PER_BASE;
		storage_type v = _bases[slot];
		v &= ~(BASE_MASK << offset);
		v |= static_cast<storage_type>(base) << offset;
		_bases[slot] = v;
		assert2((*this)[idx] == base);
	}

	friend class boost::serialization::access;

	// Serialize this BaseVec.
	template <class Archive>
	void save(Archive & ar, unsigned version) const
	{
		ar << _size;
		ar.save_binary(_bases, DIV_ROUND_UP(_size, BASES_PER_BYTE));
	}

	// Deserialize this BaseVec.
	template <class Archive>
	void load(Archive & ar, unsigned version)
	{
		ar >> _size;
		resize(_size);
		ar.load_binary(_bases, DIV_ROUND_UP(_size, BASES_PER_BYTE));
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()

	// Resizes this BaseVec to hold up to @size bases.
	void resize(size_type size)
	{
		_size = size;
		delete[] _bases;
		_bases = new storage_type[DIV_ROUND_UP(_size, BASES_PER_STORAGE_TYPE)];
	}

	// Initializes this BaseVec from a text string of A's, T's, C's, and
	// G's.
	void load_from_text(const std::string &s)
	{
		load_from_text(s.c_str(), s.size());
	}

	void load_from_text(const char text[], size_type len)
	{
		resize(len);
		for (size_type i = 0; i < len; i++)
			set(i, BaseUtils::ascii_to_bin(text[i]));
	}

	// Extracts the subsequence [beg, end] from this BaseVec and inserts it
	// into @dest.
	void extract_seq(const size_type beg,
			 const size_type end,
			 BaseVec & dest) const
	{
		size_type len;
		size_type i;
		if (end > beg) {
			len = end - beg + 1;
			dest.resize(len);
			for (i = 0; i < len; i++)
				dest.set(i, (*this)[beg + i]);
		} else {
			len = beg - end + 1;
			dest.resize(len);
			for (i = 0; i < len; i++)
				dest.set(i, (3 ^ (*this)[beg - i]));
		}
	}

	// Print the sequence contained in this BaseVec
	friend std::ostream & operator<<(std::ostream & os, const BaseVec & bv)
	{
		for (BaseVec::size_type i = 0; i < bv.size(); i++)
			os << BaseUtils::bin_to_ascii(bv[i]);
		return os;
	}

	BaseVec()
	{
		_size = 0;
		_bases = NULL;
	}

	// XXX It's currently expected that the BaseVec destructor does NOT free
	// the storage for the bases.
	~BaseVec() { }

	// Frees the storage allocated for this BaseVec.
	void destroy()
	{
		_size = 0;
		delete[] _bases;
		_bases = NULL;
	}
};
