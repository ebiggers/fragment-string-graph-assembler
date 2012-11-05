#pragma once

#include "BaseUtils.h"
#include "util.h"

#include <boost/serialization/split_member.hpp>
#include <string>
#include <ostream>
#include <assert.h>

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

	size_type size() const
	{
		return _size;
	}

	unsigned char operator[](unsigned idx) const
	{
		assert2(idx < _size);
		size_type slot = idx / BASES_PER_STORAGE_TYPE;
		size_type offset = (idx % BASES_PER_STORAGE_TYPE) * BITS_PER_BASE;
		return (_bases[slot] >> offset) & BASE_MASK;
	}

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

	template <class Archive>
	void save(Archive & ar, unsigned version) const
	{
		ar << _size;
		ar.save_binary(_bases, DIV_ROUND_UP(_size, BASES_PER_BYTE));
	}

	template <class Archive>
	void load(Archive & ar, unsigned version)
	{
		ar >> _size;
		resize(_size);
		ar.load_binary(_bases, DIV_ROUND_UP(_size, BASES_PER_BYTE));
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()

	void resize(size_type size)
	{
		_size = size;
		delete[] _bases;
		_bases = new storage_type[DIV_ROUND_UP(_size, BASES_PER_STORAGE_TYPE)];
	}

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

	void extract_seq(const size_type beg,
			 const size_type end,
			 BaseVec & dest) const
	{
		size_type len;
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

	friend std::ostream & operator<<(std::ostream & os, const BaseVec & bv)
	{
		//os << "BaseVec {_bases = \"";
		for (size_t i = 0; i < bv.size(); i++)
			os << BaseUtils::bin_to_ascii(bv[i]);
		//return os << "\"}";
		return os;
	}

	BaseVec()
	{
		_size = 0;
		_bases = NULL;
	}

	void destroy()
	{
		_size = 0;
		delete[] _bases;
		_bases = NULL;
	}
};
