#pragma once

#include "BaseUtils.h"
#include "util.h"

template <unsigned _K>
struct Kmer {
public:
	typedef unsigned long storage_type;
	typedef unsigned size_type;
private:
	static const size_type BITS_PER_BASE = 2;
	static const size_type BASES_PER_BYTE = 8 / BITS_PER_BASE;
	static const size_type BASES_PER_STORAGE_TYPE =
				sizeof(storage_type) * 8 / BITS_PER_BASE;
	static const size_type LAST_BASE_SHIFT =
				sizeof(storage_type) * 8 - BITS_PER_BASE;
	static const storage_type BASE_MASK = (1 << BITS_PER_BASE) - 1;
	static const size_type NUM_STORAGES = DIV_ROUND_UP(_K, BASES_PER_STORAGE_TYPE);

	static const storage_type LAST_STORAGE_MASK =
		static_cast<storage_type>(1) <<
			(BITS_PER_BASE *
			 	(BASES_PER_STORAGE_TYPE - (_K % BASES_PER_STORAGE_TYPE)));

	storage_type _bases[NUM_STORAGES];
public:
	static const size_type K = _K;

	Kmer() { 
		for (size_type i = 0; i < NUM_STORAGES; i++)
			_bases[i] = 0;
	}
	~Kmer() { }

	void push_back(unsigned char base)
	{
		for (size_type i = 0; i < NUM_STORAGES - 1; i++) {
			_bases[i] = (_bases[i] >> BITS_PER_BASE)
				    | _bases[i + 1] << LAST_BASE_SHIFT;
		}
		_bases[NUM_STORAGES - 1] = (_bases[NUM_STORAGES - 1] >> BITS_PER_BASE) |
					   (static_cast<storage_type>(base) << LAST_BASE_SHIFT);
	}

	void push_front(unsigned char base)
	{
		_bases[NUM_STORAGES - 1] <<= BITS_PER_BASE;
		if (NUM_STORAGES > 1) {
			size_type i = NUM_STORAGES - 1;
			do {
				storage_type last_base = _bases[i] >> LAST_BASE_SHIFT;
				_bases[i + 1] |= last_base;
				_bases[i] <<= BITS_PER_BASE;
			}  while (i-- != 0);
		}
		_bases[0] |= base;
		_bases[NUM_STORAGES - 1] |= LAST_STORAGE_MASK;
	}

	void push_rc_front(unsigned char base)
	{
		push_front(base ^ 3);
	}

	void complement()
	{
		for (size_type i = 0; i < NUM_STORAGES; i++)
			_bases[i] ^= std::numeric_limits<size_type>::max();
	}

	void reverse()
	{
		unimplemented();
	}

	void reverse_complement()
	{
		reverse();
		complement();
	}

	unsigned char operator[](unsigned idx) const
	{
		unsigned slot = idx / BASES_PER_STORAGE_TYPE;
		unsigned offset = (idx % BASES_PER_STORAGE_TYPE) * BITS_PER_BASE;
		return static_cast<unsigned char>((_bases[slot] >> offset) & BASE_MASK);
	}

	friend bool operator==(const Kmer<_K> & kmer_1, const Kmer<_K> & kmer_2)
	{
		for (size_type i = 0; i < _K; i++)
			if (kmer_1[i] != kmer_2[i])
				return false;
		return true;
	}

	friend bool operator<(const Kmer<_K> & kmer_1, const Kmer<_K> & kmer_2)
	{
		for (size_type i = 0; i < _K; i++)
			if (kmer_1[i] < kmer_2[i])
				return true;
		return false;
	}

	friend const Kmer<_K> &
	canonical_kmer(const Kmer<_K> & kmer_1, const Kmer<_K> & kmer_2)
	{
		return (kmer_1 < kmer_2) ? kmer_1 : kmer_2;
	}

	size_t hash() const
	{
		size_t h = 0;
		for (size_type i = 0; i < NUM_STORAGES; i++)
			h ^= _bases[i];
	}
	friend struct std::hash<Kmer<_K> >;
};
namespace std {
	template <unsigned _K>
	struct hash<Kmer<_K> > {
		size_t operator()(const Kmer<_K> & kmer) const
		{
			return kmer.hash();
		}
	};
}
