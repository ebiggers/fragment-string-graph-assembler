/*
 * parse-reads.cpp
 *
 * FASTA and FASTQ file parsing; interface to call a function on all the reads
 * in a set of FASTA or FASTQ files.
 */

#include <omp.h>

#include "assert.h"
#include "parse-reads.h"
#include "SharedQueue.h"
#include "Read.h"
#include "util.h"
#include "InputStream.h"
#include "OutputStream.h"

#define DEBUG(format, ...)

typedef bool (*next_read_t)(Read &r, InputStream &in, unsigned long &line_num);

static bool next_fasta_read(Read &r, InputStream &in, unsigned long &line_num)
{
	ssize_t len;

	// Read and trim the tag
	len = in.getline(&r.tag, &r.tag_buf_sz);
	if (len == -1)
		return false;
	line_num++;
	if (r.tag[0] != '>') {
		fatal_error("FASTA syntax error at line %lu of \"%s\": "
			    "expected '>' token", line_num, in.s_filename);
	}
	if (r.tag[len - 1] == '\n')
		r.tag[--len] = '\0';
	r.tag_len = len;

	// Allocate initial buffer for sequence if not already allocated
	if (!r.seq) {
		r.seq_buf_sz = 128;
		r.seq = (char*)xmalloc(r.seq_buf_sz);
	}
	len = 0;
	while (1) { // Read next FASTA sequence
		int c = in.getc();
		if (c == -1) // EOF
			goto seq_done;
		in.ungetc(c);
		if (c == '>') // Tag for next sequence
			goto seq_done;
		ssize_t line_offset = len;
		while (1) { // Read next sequence line
			char *ret = in.fgets(r.seq + len, r.seq_buf_sz - len);
			if (ret == NULL) // EOF
				goto seq_done;
			len += strlen(r.seq + len);
			if (len == (ssize_t)r.seq_buf_sz - 1) {
				// All space remaining was used.

				if (r.seq[r.seq_buf_sz - 2] == '\n')
					break; // End of line

				// Not end of line; more space needed.
				r.seq_buf_sz *= 2;
				r.seq = (char*)xrealloc(r.seq, r.seq_buf_sz);
			} else {
				// Not all the remaining space was used, so this
				// is end-of-line.
				break;
			}
		}			
		// Trim whitespace at end of line.
		len = line_offset + trim(r.seq + line_offset, len - line_offset);
		line_num++;
	}
seq_done:
	if (len == 0)
		fatal_error("FASTA syntax error at line %lu of \"%s\": "
			    "expected nonzero sequence length!", line_num,
			    in.s_filename);
	r.seq_len = len;
	//printf("%s %zu\n", r.tag, r.seq_len);
	return true;
}

static bool next_fastq_read(Read &r, InputStream &in, unsigned long &line_num)
{
	ssize_t len;
	len = in.getline(&r.tag, &r.tag_buf_sz);
	line_num++;
	if (len == -1)
		return false;
	if (r.tag[0] != '@') {
		fatal_error("FASTQ syntax error at line %lu of \"%s\": "
			    "expected '@' token", line_num, in.s_filename);
	}
	if (r.tag[len - 1] == '\n')
		r.tag[--len] = '\0';
	r.tag_len = len;

	len = in.getline(&r.seq, &r.seq_buf_sz);
	line_num++;
	if (len == -1)
		goto unexpected_eof;
	if (r.seq[len - 1] == '\n')
		r.seq[--len] = '\0';
	r.seq_len = len;
	
	len = in.getline(&r.quals, &r.quals_buf_sz);
	line_num++;
	if (len == -1)
		goto unexpected_eof;
	if (r.quals[0] != '+')
		fatal_error("FASTQ syntax error at line %lu of \"%s\": "
			    "expected '+' token", line_num, in.s_filename);

	len = in.getline(&r.quals, &r.quals_buf_sz);
	line_num++;
	if (len == -1)
		goto unexpected_eof;
	
	if (r.quals[len - 1] == '\n')
		r.quals[--len] = '\0';
	if (len != r.seq_len)
		fatal_error("FASTQ syntax error at line %lu of \"%s\": "
			    "quality string is not the same length as sequence",
			    line_num, in.s_filename);

	return true;
unexpected_eof:
	fatal_error_with_errno("Unexpected EOF in FASTQ file \"%s\"",
			       in.s_filename);
}

/* Read reads from the FASTA or FASTQ file @reads_file and add them to the
 * @ready_reads_queue queue. */
static void produce_reads(SharedQueue<ReadSet*> &free_reads_queue,
			  SharedQueue<ReadSet*> &ready_reads_queue,
			  InputStream &in,
			  unsigned &num_files_remaining, 
			  unsigned num_consumer_threads)
{
	/* Assume that if the reads file begins with the character '@' it is a
	 * FASTQ file, and if it begins with the character '>', it is a FASTA
	 * file. */
	int c = in.getc();
	in.ungetc(c);
	next_read_t next_read;
	if (c == '>')
		next_read = next_fasta_read;
	else if (c == '@')
		next_read = next_fastq_read;
	else
		fatal_error("\"%s\" does not seem to be a FASTA or FASTQ file!",
			    in.s_filename);
	
	unsigned long line_num = 0;
	unsigned long num_reads_processed = 0;
	ReadSet *s;
	unsigned i;
	unsigned n;
	while (1) {
		s = free_reads_queue.get();
		for (i = 0; i < READS_PER_SET; i++) {
			if (!next_read(s->reads[i], in, line_num))
				goto eof;
			s->reads[i].filtered = false;
			if (++num_reads_processed % 1000000 == 0) {
				info("Processed %'lu reads from file %s",
				     num_reads_processed, in.s_filename);
			}
		}
		ready_reads_queue.put(s);
	}
eof:
	// End-of-file reached!  We need to provide a partial read set.
	info("Done reading file \"%s\" (%lu reads processed)",
	     in.s_filename, num_reads_processed);
	free(s->reads[i].seq);
	s->reads[i].seq = NULL;
	ready_reads_queue.put(s);

	/* If we are the last producer to finish, poison the
	 * ready_reads_queue once for each consumer thread. */
#pragma omp critical
	n = --num_files_remaining;
	if (n == 0) {
		while (num_consumer_threads--) {
			s = free_reads_queue.get();
			s->is_last = true;
			ready_reads_queue.put(s);
		}
	}
}

struct transform_stats {
	unsigned long rd_count;
	unsigned long filtered_rd_count;
	unsigned long total_seq_len;

	transform_stats() :rd_count(0), filtered_rd_count(0), total_seq_len(0)
	{ }
};

/* Get reads from @in_queue and call @f on each one, passing it @arg, then pass
 * the reads on to @out_queue. */
static void consume_reads(SharedQueue<ReadSet*> &in_queue, 
			  SharedQueue<ReadSet*> &out_queue,
			  read_processor_t f, void *arg,
			  struct transform_stats *stats,
			  unsigned poisons_to_kill)
{

	ReadSet *s;
	while (1) {
		/* Retrieve the next read set from the queue, blocking until one
		 * is available. */
		s = in_queue.get();

		/* If this is the last set, we should return, but not before
		 * passing the set onto the next stage of consumers (or the
		 * producer's free queue, if we are the last stage of
		 * consumers). */
		if (s->is_last) {
			out_queue.put(s);
			if (--poisons_to_kill == 0)
				return;
			else
				continue;
		}

		/* Call @f on the reads in the read set. */
		for (unsigned i = 0; i < READS_PER_SET && s->reads[i].seq; i++) {
			f(s->reads[i], arg);
			if (stats) {
				stats->rd_count++;
				if (s->reads[i].filtered)
					stats->filtered_rd_count++;
				stats->total_seq_len += s->reads[i].seq_len;
			}
		}

		/* Pass the ReadSet onto the next stage of consumers, or
		 * possibly back to the producer's free_reads_queue. */
		out_queue.put(s);
	}
}

static void consume_pe_reads(SharedQueue<ReadSet*> &in_queue_1,
			     SharedQueue<ReadSet*> &in_queue_2,
			     SharedQueue<ReadSet*> &out_queue_1,
			     SharedQueue<ReadSet*> &out_queue_2,
			     Lock &unqueue_lock,
			     Lock &queue_lock,
			     pe_read_processor_t f, void *arg,
			     int lib_id, int tid)
{
	ReadSet *s1, *s2;
	while (1) {
		/* Get a set of read pairs from the reader threads.  Preserve
		 * read order! */
		unqueue_lock.lock();
		s1 = in_queue_1.get();
		s2 = in_queue_2.get();
		unqueue_lock.unlock();

		if (s1->is_last || s2->is_last) {
			/* No more reads!  Make sure that this is the case for
			 * *both* input files. */
			if (!s1->is_last || !s2->is_last) {
				fatal_error("Paired-end read files do "
					    "not have the same number "
					    "of lines");
			}
			out_queue_1.put(s1);
			out_queue_2.put(s2);
			return;
		}

		/* Call @f on each read pair in the read sets. */
		for (unsigned i = 0; i < READS_PER_SET; i++) {
			if (!s1->reads[i].seq || !s2->reads[i].seq) {
				/* One of the input files ended, so we have a
				 * partial set.  Make sure this is the case for
				 * *both* input files. */
				if (s1->reads[i].seq || s2->reads[i].seq) {
					fatal_error("Paired-end read files do "
						    "not have the same number "
						    "of lines");
				}
				break;
			}
			f(s1->reads[i], s2->reads[i], arg, lib_id, tid);
		}

		/* Hand the read set off to the writer threads. */
		queue_lock.lock();
		out_queue_1.put(s1);
		out_queue_2.put(s2);
		queue_lock.unlock();
	}
}


static void write_read(Read &r, void *arg)
{
	OutputStream &out = *(OutputStream*)arg;
	if (r.filtered)
		return;
	if (r.quals) {
		// FASTQ output
		size_t len = r.tag_len + 1 + r.seq_len + 3 + r.seq_len + 1;
		char buf[len];
		char *p = buf;
		for (ssize_t i = 0; i < r.tag_len; i++)
			*p++ = r.tag[i];
		*p++ = '\n';
		for (ssize_t i = 0; i < r.seq_len; i++)
			*p++ = r.seq[i];
		*p++ = '\n';
		*p++ = '+';
		*p++ = '\n';
		for (ssize_t i = 0; i < r.seq_len; i++)
			*p++ = r.quals[i];
		*p++ = '\n';
		out.write(buf, len);
	} else {
		// FASTA output
		size_t len = r.tag_len + 1 + r.seq_len + 1;
		char buf[len];
		char *p = buf;
		for (ssize_t i = 0; i < r.tag_len; i++)
			*p++ = r.tag[i];
		*p++ = '\n';
		for (ssize_t i = 0; i < r.seq_len; i++)
			*p++ = r.seq[i];
		*p++ = '\n';
		out.write(buf, len);
	}
}

/* Calls @f on each read in each FASTA or FASTQ file @reads_files, passing it
 * @arg as its last parameter.  Uses @num_threads threads in parallel to process
 * the reads. */
void for_all_reads(const char **reads_files, unsigned num_read_files, 
		   read_processor_t f, void *arg, unsigned num_threads)
{
	/* For each of the @num_read_files read files, create a producer thread
	 * that will read from it.
	 *
	 * The @num_threads parameter passed to the function indicates the
	 * number of consumer threads that will process the reads after being
	 * passed them by the producer thread(s).  */
	unsigned num_consumer_threads = num_threads;
	unsigned num_producer_threads = num_read_files;
	num_threads = num_consumer_threads + num_producer_threads;

	/* Make the queue size proportional to the number of producer threads.
	 * */
	size_t queue_size = num_producer_threads * QUEUE_SLOTS_PER_PRODUCER;

	/* Producers will take reads from the queue @free_reads_queue, fill them in,
	 * and add them to the queue @ready_reads_queue.  Consumers will take reads
	 * from @ready_reads_queue, process them, and return them to @free_reads_queue. 
	 *
	 * The buffers in the `class Read's, which are used to store each read's
	 * associated DNA sequence, tag, and quality scores, are re-used; the
	 * `class Read's are not destructed until the end of this function.  
	 *
	 * Initially, all the `class Read's are placed in the @free_reads_queue queue,
	 * and their buffers are all NULL (unallocated).
	 *
	 * We actually deal with sets of `class Read's (ReadSets) to improve
	 * efficiency, so we aren't passing single reads around.
	 * */
	SharedQueue<ReadSet*> free_reads_queue(queue_size);
	SharedQueue<ReadSet*> ready_reads_queue(queue_size);
	for (size_t i = 0; i < queue_size; i++)
		free_reads_queue.put(new ReadSet);

	info("Starting %u reader threads to read the following files:",
				num_producer_threads);
	for (unsigned i = 0; i < num_read_files; i++)
		info("\t\"%s\"", reads_files[i]);

	info("Starting %u consumer threads to process the reads",
				num_consumer_threads);

#pragma omp parallel num_threads(num_threads)
	{
		int thread_num = omp_get_thread_num();
		if ((unsigned)omp_get_num_threads() != num_threads)
			fatal_error("Not enough threads (is OpenMP disabled?)");

		if (thread_num < (int)num_consumer_threads) {
			// Consumer threads must be numbered [0 ... num_consumer_threads - 1]
			// because it's expected to use omp_get_thread_num() to
			// get the thread number, and each thread number must be
			// less than num_consumer_threads because that's how
			// many locks were allocated in the KmerReplacementsMap.
			consume_reads(ready_reads_queue, free_reads_queue,
				      f, arg, NULL, 1);
		} else {
			InputStream in(reads_files[thread_num - num_consumer_threads]);
			produce_reads(free_reads_queue, 
				      ready_reads_queue, 
				      in,
				      num_read_files, 
				      num_consumer_threads);
		}
	}
}

/* Calls @f on each read in the FASTA or FASTQ file @reads_file, passing it @arg
 * as its last parameter.  After calling @f on a given read, write the read to
 * @output_file, possibly after being modified in-place by @f.
 *
 * @f may also set the @filtered flag on a read, which indicates that the read
 * is not to be written. */
void transform_reads(InputStream &in, OutputStream &out,
		     read_processor_t f, void *arg)
{
	size_t queue_size = QUEUE_SLOTS_PER_PRODUCER;
	SharedQueue<ReadSet*> free_input_reads_queue(queue_size);
	SharedQueue<ReadSet*> ready_input_reads_queue(queue_size);
	SharedQueue<ReadSet*> ready_output_reads_queue(queue_size);

	for (size_t i = 0; i < queue_size; i++)
		free_input_reads_queue.put(new ReadSet);
	
	unsigned num_files_remaining = 1;

	/*
	 * We create three threads to do the read transformation.
	 * 	- Producer thread to read the FASTA or FASTQ file and produce
	 * 	  reads to the ready_input_reads_queue.  Empty reads are gotten
	 * 	  from the free_input_reads_queue.
	 * 	- Transformer thread to call the transformer function on the reads.
	 * 	  It retrieves reads from the ready_input_reads_queue,
	 * 	  transforms them, and places them in the
	 * 	  ready_output_reads_queue.
	 * 	- Writer thread to write the output FASTA or FASTQ file.  It
	 * 	  retrieves reads from the ready_output_reads_queue, writes
	 * 	  them, and returns them to the free_input_reads_queue. 
	 *
	 * So overall the reads cycle through the threads and queues as follows:
	 *
	 *       free_input_reads_queue => (producer thread) => ready_input_reads_queue =>	\
	 *       	(transformer thread) => ready_output_reads_queue =>			\
	 *       		(writer thread) => free_input_reads_queue
	 */      		
#pragma omp parallel num_threads(3)
	{
		int thread_num = omp_get_thread_num();
		if ((unsigned)omp_get_num_threads() != 3)
			fatal_error("Not enough threads (is OpenMP disabled?)");
		if (thread_num == 0) {
			produce_reads(free_input_reads_queue,
				      ready_input_reads_queue,
				      in, num_files_remaining, 1);
		} else if (thread_num == 1) {
			consume_reads(ready_input_reads_queue,
				      ready_output_reads_queue, f, arg,
				      NULL, 1);
		} else {
			consume_reads(ready_output_reads_queue,
				      free_input_reads_queue,
				      write_read, &out, NULL, 1);
		}
	}
	info("Done writing file \"%s\"", out.s_filename);
}

/*
 * Processing a single library for the transform_pe_reads_files() function.
 *
 * The idea here is to have a thread to read each of the two FASTA or FASTQ
 * files, a thread to write each of the two output FASTA or FASTQ files, and an
 * arbitrary number of transformer threads that repeatedly get a read pair from
 * the reader threads, do some sort of read modification or filtering, then hand
 * the read pair off to the writer threads.
 */
static void transform_pe_reads(InputStream &in_1, InputStream &in_2,
			       OutputStream &out_1, OutputStream &out_2,
			       pe_read_processor_t f, void *arg,
			       int num_transformer_threads,
			       struct transform_stats *stats_1,
			       struct transform_stats *stats_2,
			       int lib_id)
{
	/*  ------>-------                              ---------<---------
	 *  |            |                              |                 |
	 *  |            V                              V                 |
	 *  |  free_input_reads_queue_1     free_input_reads_queue_2      |
	 *  |            |                              |                 |
	 *  |            V                              V                 |
	 *  |  %%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%% |
	 *  ^  %       THREAD            %    %         THREAD          % ^
	 *  |  %    Read reader 1        %    %     Read reader 2       % |
	 *  |  % reading from file @in_1 %    % reading from file @in_2 % |
	 *  |  %%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%% |
	 *  |            |                              |                 |
	 *  |            V                              V                 |
	 *  |  ready_input_reads_queue_1    ready_input_reads_queue_2     |
	 *  |            |                              |                 |
	 *  |            ---------------    -------------                 |
	 *  |     ...                  |    |                 ...         |
	 *  ^                          V    V                             ^
	 *  | %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%  |
	 *  | %    THREAD     % %     THREAD    %      %    THREAD     %  | 
	 *  | % transformer 1 % % transformer 2 % ...  % transformer n %  |
	 *  | %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%  |
	 *  |                          |    |                             |
	 *  |     ...    ---------------    -------------       ...       |
	 *  |            |                              |                 |
	 *  |            V                              V                 |
	 *  ^  ready_output_reads_queue_1    ready_output_reads_queue_2   ^
	 *  |            |                              |                 |
	 *  |            V                              V                 |
	 *  |  %%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%  |
	 *  |  %       THREAD           %     %       THREAD           %  |
	 *  |  %     Read writer 1      %     %   Read writer 2        %  |
	 *  |  % writing to file @out_1 %     % writing to file @out_2 %  |
	 *  |  %%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%  |
	 *  |            |                              |                 |
	 *  ----<-----<---                              ------>-------->---
	 *
	 */

	/* Caution: it's important for ordering of the reads to be maintained,
	 * so the shared queues must preserve the ordering property, and the
	 * un-queueing and queueing of the two reads in a read pair must be
	 * protected by a lock so that they stay in order.
	 *
	 * But because these locks must exist, get() and put() are each only
	 * ever called by 1 thread at a time on all these shared queues, so it
	 * is indeed safe to use the lock-free version of the SharedQueue (see
	 * SharedQueue.h), even though it cannot guarantee the ordering property
	 * in the event of concurrent get()s or put()s. */
	size_t queue_size = QUEUE_SLOTS_PER_PRODUCER;
	SharedQueue<ReadSet*> free_input_reads_queue_1(queue_size);
	SharedQueue<ReadSet*> free_input_reads_queue_2(queue_size);
	SharedQueue<ReadSet*> ready_input_reads_queue_1(queue_size);
	SharedQueue<ReadSet*> ready_input_reads_queue_2(queue_size);
	SharedQueue<ReadSet*> ready_output_reads_queue_1(queue_size);
	SharedQueue<ReadSet*> ready_output_reads_queue_2(queue_size);

	for (size_t i = 0; i < queue_size; i++) {
		free_input_reads_queue_1.put(new ReadSet);
		free_input_reads_queue_2.put(new ReadSet);
	}

	Lock unqueue_lock;
	Lock queue_lock;
	int num_threads = 4 + num_transformer_threads;
#pragma omp parallel num_threads(num_threads)
	{
		int thread_num = omp_get_thread_num();
		if (omp_get_num_threads() != num_threads)
			fatal_error("Not enough threads (is OpenMP disabled?)");
		if (thread_num < num_transformer_threads) {
			// Transformer thread
			consume_pe_reads(ready_input_reads_queue_1,
					 ready_input_reads_queue_2,
					 ready_output_reads_queue_1,
					 ready_output_reads_queue_2,
					 unqueue_lock, queue_lock, f, arg,
					 lib_id, thread_num);
		} else if (thread_num == num_transformer_threads) {
			// PE read reader thread 1
			unsigned num_files_remaining = 1;
			produce_reads(free_input_reads_queue_1,
				      ready_input_reads_queue_1,
				      in_1, num_files_remaining,
				      num_transformer_threads);
		} else if (thread_num == num_transformer_threads + 1) {
			// PE read reader thread 2
			unsigned num_files_remaining = 1;
			produce_reads(free_input_reads_queue_2,
				      ready_input_reads_queue_2,
				      in_2, num_files_remaining,
				      num_transformer_threads);
		} else if (thread_num == num_transformer_threads + 2) {
			// PE read writer thread 1
			consume_reads(ready_output_reads_queue_1,
				      free_input_reads_queue_1,
				      write_read, &out_1, stats_1,
				      num_transformer_threads);
		} else if (thread_num == num_transformer_threads + 3) {
			// PE read writer thread 2
			consume_reads(ready_output_reads_queue_2,
				      free_input_reads_queue_2,
				      write_read, &out_2, stats_2,
				      num_transformer_threads);
		} else {
			assert(0);
		}
	}
}

static void print_transform_info(const char *in_filename,
				 const char *out_filename,
				 const struct transform_stats &stats)
{
	info("Done transforming \"%s\" => \"%s\"", in_filename, out_filename);
	info("  Number of reads processed: %'lu", stats.rd_count);
	if (stats.rd_count) {
		info("  Number of reads filtered: %'lu (%.2f%%)",
		     stats.filtered_rd_count,
		     double(stats.filtered_rd_count * 100) / stats.rd_count);
		info("  Average read length: %lu bp",
		     stats.total_seq_len / stats.rd_count);
	}
}

/*
 * Calls a function on each mate pair in a list of paired-end read libraries,
 * then writes the mate pair to new FASTA or FASTQ files, unless it was marked
 * as filtered.
 *
 * This is done in parallel across libraries and within each library.  Although
 * all the reads will be paired up properly, there is no guarantee as to the
 * ordering of the mate pairs.
 *
 *				  Parameters:
 *
 * @reads_files:    The filenames of FASTQ or FASTA files containing the
 *                  paired-end-reads, where each consecutive pair of filenames
 *                  gives read 1 and read 2 of a library.
 *
 * @num_reads_files:  Number of filenames provided (must be even)
 *
 * @added_suffix:  A suffix to add to the filenames when they're transformed.
 * 		   The suffix is added *before* the fq, fastq, or fa suffix, if
 * 		   present, and the transformed files are placed in the working
 * 		   directory.
 *
 * @f:	A function that will be called on each mate pair in the libraries.
 *      f may modify each read in-place, and it also may set the @filtered flag
 *      on a read to cause that read to not be written.  (But if you filter one
 *      read in a mate pair, you really should filter both.)
 *
 * @args:	Additional argument to pass @f.
 *
 * @num_threads:  Number of transformer threads per paired-end library.
 */
void transform_pe_reads_files(const char          **reads_files,
			      int		    num_reads_files,
			      const char	   *added_suffix,
			      pe_read_processor_t   f,
			      void		   *arg,
			      int		    num_threads)
{
	if (num_reads_files % 2 != 0)
		fatal_error("Expected even number of paired-end reads files");

	omp_set_nested(true);

	int num_libs = num_reads_files / 2;
	struct transform_stats stats[num_reads_files];
	Lock print_lock;
	#pragma omp parallel for num_threads(num_libs)
	for (int i = 0; i < num_libs; i++) {
		char *out_reads_1, *out_reads_2;
		char *log_1, *log_2;
		get_output_filenames(reads_files[i * 2], added_suffix,
				     &out_reads_1, &log_1);
		get_output_filenames(reads_files[i * 2 + 1], added_suffix,
				     &out_reads_2, &log_2);
		InputStream in_1(reads_files[i * 2]);
		InputStream in_2(reads_files[i * 2 + 1]);
		OutputStream out_1(out_reads_1);
		OutputStream out_2(out_reads_2);

		print_lock.lock();
		info("Transforming paired-end library:");
		info(" (\"%s\", \"%s\")", reads_files[i * 2], reads_files[i * 2 + 1]);
		info("  => (\"%s\", \"%s\")", out_reads_1, out_reads_2);
		info(" ");
		print_lock.unlock();

		transform_pe_reads(in_1, in_2, out_1, out_2, f, arg,
				   num_threads, &stats[i * 2],
				   &stats[i * 2 + 1], i);

		print_lock.lock();
		print_transform_info(reads_files[i * 2],
				     out_reads_1, stats[i * 2]);

		print_transform_info(reads_files[i * 2 + 1],
				     out_reads_2, stats[i * 2 + 1]);
		print_lock.unlock();

		free(out_reads_1);
		free(out_reads_2);
		free(log_1);
		free(log_2);
	}
}
