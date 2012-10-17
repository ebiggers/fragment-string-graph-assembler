/*
 * dna-util.h
 *
 * Functions for manipulating DNA bases/characters and k-mers.
 */

#ifndef _DNA_UTIL_H
#define _DNA_UTIL_H

#include <inttypes.h>

extern uint64_t reverse_complement(uint64_t kmer, unsigned kmer_len);
extern const char *kmer_str(uint64_t kmer, unsigned kmer_len);
extern void kmer_str_r(uint64_t kmer, unsigned kmer_len, char buf[]);
extern const unsigned char ascii_to_bin_tab[256];
extern const char bin_to_ascii_tab[4];
extern const char ascii_complement_tab[256];


/* Returns the canonical form of two k-mers. */
static inline uint64_t canonical_form(uint64_t kmer_1, uint64_t kmer_2)
{
	return (kmer_1 < kmer_2) ? kmer_1 : kmer_2;
}

/* Returns the canonical form of a k-mer. */
static inline uint64_t to_canonical_form(uint64_t kmer, unsigned kmer_len)
{
	return canonical_form(kmer, reverse_complement(kmer, kmer_len));
}

/* Returns the complement of a binary base. */
static inline unsigned char complement(unsigned char base)
{
	return base ^ 3;
}

static inline char ascii_complement(char base)
{
	return ascii_complement_tab[(unsigned char)base];
}

/* Returns true iff the kmer is in canonical form. */
static inline bool is_canonical(uint64_t kmer, unsigned kmer_len)
{
	return to_canonical_form(kmer, kmer_len) == kmer;
}

/* Translates a nucleobase from ASCII to binary. */
static inline unsigned char base_ascii_to_bin(unsigned char c)
{
	return ascii_to_bin_tab[c];
}

/* Translates a nucleobase from ASCII to binary (signed char version) */
static inline unsigned char base_ascii_to_bin(char c)
{
	return ascii_to_bin_tab[(unsigned char)c];
}

/* Translates a nucleobase from binary to ASCII */
static inline char base_bin_to_ascii(unsigned char base)
{
	return bin_to_ascii_tab[base];
}

/* Returns true if an `unsigned char' number is a valid nucleobase in the binary
 * representation. */
static inline bool bin_base_is_dna(unsigned char base)
{
	return (base < 4);
}


/* Adds a binary nucleobase to the end of a sequence of nucleobases stored
 * together in a 64-bit number. */
static inline uint64_t mer_add_base(uint64_t mer, unsigned char base)
{
	return (mer << 2) | base;
}

/* Adds a binary nucleobase to the end of a sequence of nucleobases stored
 * together in a 64-bit number that is @kmer_len bases long, while removing the
 * first base. */
static inline uint64_t kmer_shift_base(uint64_t kmer, unsigned char base,
				       unsigned kmer_len)
{
	kmer = mer_add_base(kmer, base);
	return (kmer & (((uint64_t)1 << (kmer_len << 1)) - 1));
}

/* Adds a binary nucleobase to the beginning of a sequence of nucleobases stored
 * together in a 64-bit number that is @kmer_len bases long.  The last base is
 * removed to make room. */
static inline uint64_t kmer_reverse_shift_base(uint64_t kmer, 
					       unsigned char base, 
					       unsigned kmer_len)
{
	return (kmer >> 2) | ((uint64_t)base << ((kmer_len - 1) << 1));
}

/* Returns the middle nucleobase of a kmer. */
static inline unsigned char kmer_middle_base(uint64_t kmer, unsigned kmer_len)
{
	return (kmer >> (kmer_len - 1)) & 3;
}

/* Returns true iff the middle nucleobase of the specified kmer is an 'A'. */
static inline bool kmer_middle_base_is_adenine(uint64_t kmer, unsigned kmer_len)
{
	return (kmer & ((uint64_t)3 << (kmer_len - 1))) == 0;
}

/* Replaces the middle nucleobase of the specified kmer with the binary
 * nucleobase @base. */
static inline uint64_t kmer_replace_middle_base(uint64_t kmer, 
						unsigned char base,
						unsigned kmer_len)
{
	unsigned shift = kmer_len - 1;
	uint64_t mask = ~((uint64_t)3 << shift);
	return (kmer & mask) | ((uint64_t)base << shift);
}

/* Returns true iff the k-mers @kmer_1 and @kmer_2 of length @kmer_len are
 * equal, disregarding the middle base. */
static inline bool kmers_equal_except_middle_base(uint64_t kmer_1, 
						  uint64_t kmer_2,
						  unsigned kmer_len)
{
	unsigned shift = kmer_len - 1;
	uint64_t mask = ~((uint64_t)3 << shift);
	return (kmer_1 & mask) == (kmer_2 & mask);
}

static inline bool kmers_equal_except_middle_region(uint64_t kmer_1, 
						    uint64_t kmer_2,
						    unsigned kmer_len)
{
	unsigned shift1 = kmer_len;
	unsigned shift3 = kmer_len * 3;
	uint64_t mask1 = ((uint64_t)1 << shift1) - 1;
	uint64_t mask3 = ((uint64_t)1 << shift3) - 1;
	uint64_t mask23 = mask3 ^ mask1;
	return (kmer_1 & ~mask23) == (kmer_2 & ~mask23);
}

/* Returns the @idx-th base of the k-mer @kmer of length @kmer_len, where the
 * bases are numbered starting at 0. */
static inline unsigned char kmer_get_base(uint64_t kmer, unsigned idx, 
					  unsigned kmer_len)
{
	unsigned shift = (kmer_len - idx - 1) << 1;
	return (kmer >> shift) & 3;
}

/* Replaces the @idx-th base of the k-mer @kmer of length @kmer_len with the new
 * binary base @base, where the bases are numbered starting at 0. */
static inline uint64_t kmer_replace_base(uint64_t kmer, unsigned char base,
					 unsigned idx, unsigned kmer_len)
{
	unsigned shift = (kmer_len - idx - 1) << 1;
	uint64_t mask = ~((uint64_t)3 << shift);
	return (kmer & mask) | ((uint64_t)base << shift);
}

#endif
