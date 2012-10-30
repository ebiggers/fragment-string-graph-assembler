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
	static const size_type BITS_PER_BYTE = 8;
	static const size_type BASES_PER_BYTE = BITS_PER_BYTE / BITS_PER_BASE;
	static const storage_type BASE_MASK =
			(static_cast<storage_type>(1) << BITS_PER_BASE) - 1;

	static const size_type BASES_PER_STORAGE_TYPE =
				sizeof(storage_type) * BASES_PER_BYTE;

	static const size_type NUM_STORAGES =
				DIV_ROUND_UP(_K, BASES_PER_STORAGE_TYPE);

	static const size_type BASES_IN_LAST_STORAGE =
				MODULO_NONZERO(_K, BASES_PER_STORAGE_TYPE);

	static const size_type STORAGE_LAST_BASE_SHIFT =
				sizeof(storage_type) * BITS_PER_BYTE - BITS_PER_BASE;

	static const size_type LAST_BASE_SHIFT =
				(BASES_IN_LAST_STORAGE - 1) * BITS_PER_BASE;

	static const storage_type LAST_STORAGE_MASK =
		(static_cast<storage_type>(1) <<
		 	(BASES_IN_LAST_STORAGE * BITS_PER_BASE)) - 1;

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
		if (NUM_STORAGES >= 2) {
			size_type i = 0;
			do {
				_bases[i] = (_bases[i] >> BITS_PER_BASE)
					    | _bases[i + 1] << STORAGE_LAST_BASE_SHIFT;
			} while (i++ != NUM_STORAGES - 2);
		}
		_bases[NUM_STORAGES - 1] = (_bases[NUM_STORAGES - 1] >> BITS_PER_BASE) |
					   (static_cast<storage_type>(base) << LAST_BASE_SHIFT);
	}

	void push_front(unsigned char base)
	{
		_bases[NUM_STORAGES - 1] <<= BITS_PER_BASE;
		if (NUM_STORAGES >= 2) {
			size_type i = NUM_STORAGES - 2;
			do {
				storage_type last_base = _bases[i] >> STORAGE_LAST_BASE_SHIFT;
				_bases[i + 1] |= last_base;
				_bases[i] <<= BITS_PER_BASE;
			}  while (i-- != 0);
		}
		_bases[0] |= base;
		_bases[NUM_STORAGES - 1] &= LAST_STORAGE_MASK;
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
		for (size_type i = 0; i < Kmer<_K>::NUM_STORAGES; i++)
			if (kmer_1._bases[i] != kmer_2._bases[i])
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
		//static const uint64_t GOLDEN_RATIO_PRIME_64 = 0x9e37fffffffc0001UL;
		size_t h = 14695981039346656037ul;
		for (size_type i = 0; i < NUM_STORAGES; i++)
			h = 1099511628211ul * (h ^ _bases[i]);
		return h;
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
