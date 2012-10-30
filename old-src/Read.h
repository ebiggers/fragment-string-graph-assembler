/*
 * Read.h
 */

#ifndef _READ_H
#define _READ_H

#include <stddef.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include "util.h"

#define READS_PER_SET 32
#define QUEUE_SLOTS_PER_PRODUCER 32

/* Represents a DNA read from a FASTA or FASTQ file. */
class Read {
public:
	char	*tag;
	char	*seq;
	char	*quals;
	ssize_t  tag_len;
	ssize_t	 seq_len;
	size_t	 tag_buf_sz;
	size_t	 seq_buf_sz;
	size_t	 quals_buf_sz;
	bool     filtered;

	Read() : tag(NULL), 
		 seq(NULL), 
		 quals(NULL), 
		 seq_len(0), 
		 tag_buf_sz(0),
		 seq_buf_sz(0), 
		 quals_buf_sz(0) 
		{ }

	Read(const Read &other)
	{
		tag   = (char*)xmalloc(other.tag_len);
		seq   = (char*)xmalloc(other.seq_len);
		quals = (char*)xmalloc(other.seq_len);
		memcpy(tag, other.tag, other.tag_len);
		memcpy(seq, other.seq, other.seq_len);
		memcpy(quals, other.quals, other.seq_len);
		tag_len      = other.tag_len;
		seq_len      = other.seq_len;
		tag_buf_sz   = other.tag_len;
		seq_buf_sz   = other.seq_len;
		quals_buf_sz = other.seq_len;
		filtered     = false;
	}

	~Read() { 
		free(tag); 
		free(seq); 
		free(quals); 
	}
};


class ReadSet {
public:
	Read reads[READS_PER_SET];
	bool is_last;

	ReadSet() : is_last(false){ }
};

#endif
