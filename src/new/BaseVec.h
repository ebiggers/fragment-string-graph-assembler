#pragma once

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/binary_object.hpp>
#include "BaseUtils.h"

class BaseVec {
public:
	typedef unsigned char storage_type;
	typedef unsigned size_type;
private:
	size_type _size;
	storage_type *_bases;
	static const size_type BITS_PER_BASE = 2;
	static const size_type BASES_PER_STORAGE_TYPE = sizeof(storage_type) * 8 / BITS_PER_BASE;
	static const size_type BASE_MASK = (1 << BITS_PER_BASE) - 1;
public:

	size_type size() const {
		return _size;
	}

	unsigned char operator[](size_type idx) const {
		size_type slot = idx / BASES_PER_STORAGE_TYPE;
		size_type offset = (idx % BASES_PER_STORAGE_TYPE) * 2;
		return (_bases[slot] >> offset) & BASE_MASK;
	}

	void set(size_type idx, storage_type base) {
		size_type slot = idx / BASES_PER_STORAGE_TYPE;
		size_type offset = (idx % BASES_PER_STORAGE_TYPE) * 2;
		storage_type v = _bases[slot];
		v &= ~(BASE_MASK << offset);
		v |= base << offset;
		_bases[slot] = v;
	}

	template <class Archive>
	void save(Archive & ar, unsigned version) {
		ar << _size;
		ar << boost::serialization::binary_object(_bases, _size);
	}

	template <class Archive>
	void load(Archive & ar, unsigned version) {
		ar >> _size;
		resize(_size);
		ar >> boost::serialization::binary_object(_bases, _size);
	}

	void resize(size_type size) {
		_size = size;
		delete[] _bases;
		_bases = new storage_type[(_size + BASES_PER_STORAGE_TYPE - 1) /
								  BASES_PER_STORAGE_TYPE];
	}

	void load_from_text(const char text[], size_type len) {
		resize(len);
		for (size_type i = 0; i < len; i++) {
			set(i, BaseUtils::ascii_to_bin(text[i]));
		}
	}

	BaseVec() {
		_size = 0;
		_bases = NULL;
	}

	~BaseVec() {
		_size = 0;
		delete[] _bases;
	}
};
