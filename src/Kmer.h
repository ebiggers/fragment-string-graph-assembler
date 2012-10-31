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

	static const size_type BASES_PER_STORAGE =
				sizeof(storage_type) * BASES_PER_BYTE;

	static const size_type NUM_STORAGES =
				DIV_ROUND_UP(_K, BASES_PER_STORAGE);

	static const size_type BASES_IN_PARTIAL_STORAGE =
				MODULO_NONZERO(_K, BASES_PER_STORAGE);

	static const storage_type PARTIAL_STORAGE_MASK =
		(BASES_IN_PARTIAL_STORAGE == BASES_PER_STORAGE) ?
			std::numeric_limits<storage_type>::max() :
			(static_cast<storage_type>(1) <<
				(BASES_IN_PARTIAL_STORAGE * BITS_PER_BASE)) - 1;

	static const size_type FIRST_BASE_SHIFT =
				(BASES_PER_STORAGE - 1) * BITS_PER_BASE;

	static const size_type PARTIAL_STORAGE_FIRST_BASE_SHIFT =
				(BASES_IN_PARTIAL_STORAGE - 1) * BITS_PER_BASE;


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
		_bases[0] = (_bases[0] << BITS_PER_BASE) & PARTIAL_STORAGE_MASK;
		if (NUM_STORAGES >= 2) {
			size_type i = 0;
			do {
				_bases[i] = _bases[i] | (_bases[i + 1] >> FIRST_BASE_SHIFT);
				_bases[i + 1] <<= BITS_PER_BASE;
			} while (++i != NUM_STORAGES - 1);
		}
		_bases[NUM_STORAGES - 1] |= base;
	}

	void push_front(unsigned char base)
	{
		if (NUM_STORAGES >= 2) {
			size_type i = NUM_STORAGES - 1;
			do {
				_bases[i] = (_bases[i] >> BITS_PER_BASE)
						| _bases[i - 1] << FIRST_BASE_SHIFT;
			}  while (--i != 0);
		}
		_bases[0] = (_bases[0] >> BITS_PER_BASE) |
			    (static_cast<storage_type>(base) << PARTIAL_STORAGE_FIRST_BASE_SHIFT);
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

	unsigned char operator[](const unsigned idx) const
	{
		unsigned slot = (idx + BASES_PER_STORAGE - BASES_IN_PARTIAL_STORAGE) /
				       BASES_PER_STORAGE;
		unsigned shift = ((BASES_PER_STORAGE - 1) -
				 ((idx + BASES_PER_STORAGE - BASES_IN_PARTIAL_STORAGE) %
				 BASES_PER_STORAGE)) * BITS_PER_BASE;
		return static_cast<unsigned char>((_bases[slot] >> shift) & BASE_MASK);
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
		for (size_type i = 0; i < NUM_STORAGES; i++)
			if (kmer_1._bases[i] < kmer_2._bases[i])
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
