#pragma once

#include "BaseUtils.h"
#include "util.h"

#include <boost/serialization/split_member.hpp>
#include <string>

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
				sizeof(storage_type) * 8 / BITS_PER_BASE;
	static const storage_type BASE_MASK = (1 << BITS_PER_BASE) - 1;

public:

	size_type size() const
	{
		return _size;
	}

	unsigned char operator[](unsigned idx) const
	{
		unsigned slot = idx / BASES_PER_STORAGE_TYPE;
		unsigned offset = (idx % BASES_PER_STORAGE_TYPE) * BITS_PER_BASE;
		return static_cast<unsigned char>((_bases[slot] >> offset) & BASE_MASK);
	}

	void set(size_type idx, unsigned char base)
	{
		size_type slot = idx / BASES_PER_STORAGE_TYPE;
		size_type offset = (idx % BASES_PER_STORAGE_TYPE) *
				    BITS_PER_BASE;
		storage_type v = _bases[slot];
		v &= ~(BASE_MASK << offset);
		v |= static_cast<storage_type>(base) << offset;
		_bases[slot] = v;
	}

	friend class boost::serialization::access;

	template <class Archive>
	void save(Archive & ar, unsigned version) const
	{
		ar << _size;
		ar.save_binary(_bases,
			       (_size + BASES_PER_BYTE - 1) / BASES_PER_BYTE);
	}

	template <class Archive>
	void load(Archive & ar, unsigned version)
	{
		ar >> _size;
		resize(_size);
		ar.load_binary(_bases,
			       (_size + BASES_PER_BYTE - 1) / BASES_PER_BYTE);
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
