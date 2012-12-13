#pragma once

#include "BaseUtils.h"
#include "util.h"

#include <iostream>
#include <BaseUtils.h>

//
// A sequence of _K bases, stored in binary format (2 bits per base).
//
// Somewhat similar to a BaseVec, but this is templatized by the length _K so it
// takes a constant amount of space determined at template instantiation time.
//
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
			~(storage_type)0 :
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

	// Push a base onto the end of the k-mer and shift all the other bases
	// left by 1 space.  The base at the front is discarded.
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

	// Push a base onto the front of the k-mer and shift all the other bases
	// right by 1 space.  The base at the end is discarded.
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

	// Changes this k-mer to the complement sequence.
	void complement()
	{
		_bases[0] ^= PARTIAL_STORAGE_MASK;
		for (size_type i = 1; i < NUM_STORAGES; i++)
			_bases[i] ^= std::numeric_limits<size_type>::max();
	}

	// Changes this k-mer to the reverse sequence.
	void reverse()
	{
		unimplemented();
	}

	// Changes this k-mer to the reverse-complement sequence.
	void reverse_complement()
	{
		reverse();
		complement();
	}

	// Return the binary base at index @idx of the k-mer.
	unsigned char operator[](const unsigned idx) const
	{
		unsigned slot = (idx + BASES_PER_STORAGE - BASES_IN_PARTIAL_STORAGE) /
				       BASES_PER_STORAGE;
		unsigned shift = ((BASES_PER_STORAGE - 1) -
				 ((idx + BASES_PER_STORAGE - BASES_IN_PARTIAL_STORAGE) %
				 BASES_PER_STORAGE)) * BITS_PER_BASE;
		return static_cast<unsigned char>((_bases[slot] >> shift) & BASE_MASK);
	}

	// Return true iff two k-mers are equal base-for-base.
	friend bool operator==(const Kmer<_K> & kmer_1, const Kmer<_K> & kmer_2)
	{
		assert2((kmer_1._bases[0] & ~PARTIAL_STORAGE_MASK) == 0);
		assert2((kmer_2._bases[0] & ~PARTIAL_STORAGE_MASK) == 0);
		for (size_type i = 0; i < NUM_STORAGES; i++)
			if (kmer_1._bases[i] != kmer_2._bases[i])
				return false;
		return true;
	}

	// Return true iff the first k-mer is lexicographically less than the
	// second k-mer.
	friend bool operator<(const Kmer<_K> & kmer_1, const Kmer<_K> & kmer_2)
	{
		assert2((kmer_1._bases[0] & ~PARTIAL_STORAGE_MASK) == 0);
		assert2((kmer_2._bases[0] & ~PARTIAL_STORAGE_MASK) == 0);
		for (size_type i = 0; i < NUM_STORAGES; i++)
			if (kmer_1._bases[i] != kmer_2._bases[i])
				return kmer_1._bases[i] < kmer_2._bases[i];
		return false;
	}

	// Returns the lexicographically lesser of two k-mers.
	friend const Kmer<_K> &
	canonical_kmer(const Kmer<_K> & kmer_1, const Kmer<_K> & kmer_2)
	{
		return (kmer_1 < kmer_2) ? kmer_1 : kmer_2;
	}

	struct hash_functor
	{
		size_t operator()(const Kmer<_K> & kmer) const
		{
			return kmer.hash();
		}
	};

	// Hashes a k-mer from its bases.
	size_t hash() const
	{
		//static const uint64_t GOLDEN_RATIO_PRIME_64 = 0x9e37fffffffc0001UL;
		assert((_bases[0] & ~PARTIAL_STORAGE_MASK) == 0);
		size_t h = 14695981039346656037ul;
		for (size_type i = 0; i < NUM_STORAGES; i++)
			h = 1099511628211ul * (h ^ _bases[i]);
		return h;
	}

	friend std::ostream & operator<<(std::ostream & os, const Kmer<_K> & kmer)
	{
		for (unsigned i = 0; i < _K; i++)
			os << BaseUtils::bin_to_ascii(kmer[i]);
	}
};
