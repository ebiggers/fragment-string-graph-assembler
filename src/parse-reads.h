/*
 * parse-reads.h
 *
 * Interface to call a function on all the reads in a set of FASTA or FASTQ
 * files.
 */

#ifndef _PARSE_READS_H
#define _PARSE_READS_H

class Read;
class InputStream;
class OutputStream;

typedef void (*read_processor_t)(Read &r, void *arg);

typedef void (*pe_read_processor_t)(Read &r1, Read &r2, void *arg,
				    int lib_id, int tid);

extern void for_all_reads(const char **reads_files, unsigned num_read_files, 
			  read_processor_t f, void *arg, unsigned num_threads);

extern void transform_reads(InputStream &in, OutputStream &out,
			    read_processor_t f, void *arg);

extern void transform_pe_reads_files(const char **reads_files,
				     int num_reads_files,
				     const char *added_suffix,
				     pe_read_processor_t f,
				     void *arg,
				     int num_threads);


#endif
