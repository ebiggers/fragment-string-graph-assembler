/*
 * util.cpp
 *
 * Miscellaneous error, warning, file, memory, and threading functions.
 */

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#include "dna-util.h"
#include "util.h"
#include "Lock.h"

static Lock output_lock;

/* Prints an error message and exits the program with failure status. */
void fatal_error(const char *msg, ...)
{
	va_list va;
	output_lock.lock();
	va_start(va, msg);
	fflush(stdout);
	fputs("ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
	output_lock.unlock();
	exit(1);
}

void fatal_error_with_errno(const char *msg, ...)
{
	va_list va;
	output_lock.lock();
	va_start(va, msg);
	fflush(stdout);
	fputs("ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	fprintf(stderr, ": %s\n", strerror(errno));
	va_end(va);
	output_lock.unlock();
	exit(1);
}

/* Prints a warning message. */
void warning(const char *msg, ...)
{
	va_list va;
	output_lock.lock();
	va_start(va, msg);
	fputs("WARNING: ", stderr);
	fprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
	output_lock.unlock();
}

const char *info_tag = NULL;

void info(const char *format, ...)
{
	va_list va;
	struct tm *tm;
	time_t t;
	static const char *month_names[] =
		{"Jan", "Feb", "Mar", "Apr", "May", "June",
		 "July", "Aug", "Sept", "Oct", "Nov", "Dec"};

	output_lock.lock();

	t = time(NULL);

	tm = localtime(&t);
	va_start(va, format);
	if (info_tag)
		fprintf(stdout, "[%s] ", info_tag);
	fprintf(stdout, "%s %d %d %2d:%02d:%02d: ",
		month_names[tm->tm_mon],
		tm->tm_mday,
		1900 + tm->tm_year,
		tm->tm_hour, tm->tm_min, tm->tm_sec);
	if (format)
		vfprintf(stdout, format, va);
	fputc('\n', stdout);
	fflush(stdout);
	va_end(va);
	output_lock.unlock();
}

/* Returns the number of processors (if it can be determined), otherwise returns
 * 1. */
unsigned get_default_num_threads()
{
	long nthreads = sysconf(_SC_NPROCESSORS_ONLN);
	if (nthreads == -1) {
		warning("Could not deteremine number of processors! Assuming 1");
		return 1;
	} else {
		return (unsigned)nthreads;
	}
}

/* Returns the reverse complement of a kmer */
uint64_t reverse_complement(uint64_t kmer, unsigned kmer_len)
{
	kmer = ((kmer >> 2)  & 0x3333333333333333ULL) | 
	    ((kmer & 0x3333333333333333ULL) << 2);
	kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FULL) | 
	    ((kmer & 0x0F0F0F0F0F0F0F0FULL) << 4);
	kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFULL) | 
	    ((kmer & 0x00FF00FF00FF00FFULL) << 8);
	kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFULL) | 
	    ((kmer & 0x0000FFFF0000FFFFULL) << 16);
	kmer = (kmer >> 32) | (kmer << 32);
	return ~kmer >> (64 - (kmer_len << 1));
}


/* Map a DNA character, upper case or lower case, to the binary code
 * { A => 0, C => 1, G => 2, C => 3, other => 4 }.
 *
 * C99 designated initializers do not exist in C++, too bad...
 */
const unsigned char ascii_to_bin_tab[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, /* A */ 0, 4, /* C */ 1, 4, 4, 4, /* G */ 2,
                4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, /* T */ 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, /* a */ 0, 4, /* c */ 1, 4, 4, 4, /* g */ 2,
                4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, /* t */ 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

/* Map a nucleotide in binary form to ASCII. */
const char bin_to_ascii_tab[4] = { 'A', 'C', 'G', 'T' };

/* Translates a kmer into the ASCII string representing it.  Warning:
 * not re-entrant */
const char *kmer_str(uint64_t kmer, unsigned kmer_len)
{
	static char buf[33];
	kmer_str_r(kmer, kmer_len, buf);
	return buf;
}

/* Re-entrant version of kmer_str().  buf should be at least (kmer_len + 1)
 * characters long. */
void kmer_str_r(uint64_t kmer, unsigned kmer_len, char buf[])
{
	char *p = buf + kmer_len;
	*p-- = '\0';
	do {
		*p = base_bin_to_ascii(kmer & 3);
		kmer >>= 2;
	} while (p-- != buf);
}

void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p) {
		fatal_error("Out of memory: tried to allocate %lu bytes",
			     size);
	}
	return p;
}

void *xzalloc(size_t size)
{
	return memset(xmalloc(size), 0, size);
}

void *xrealloc(void *ptr, size_t size)
{
	void *p = realloc(ptr, size);
	if (!p)
		fatal_error("Out of memory: tried to reallocate %lu bytes", size);
	return p;
}

void *xmemalign(size_t alignment, size_t size)
{
	void *p;
	int ret = posix_memalign(&p, alignment, size);
	if (ret != 0)
		fatal_error("Out of memory: tried to allocate %lu bytes "
				"on %zu byte boundary", size, alignment);
	return p;
}

char *xstrdup(const char *s)
{
	size_t len = strlen(s);
	char *p = (char*)xmalloc(len + 1);
	return (char*)memcpy(p, s, len + 1);
}

char *xstrdup2(const char *s1, const char *s2)
{
	size_t len1 = strlen(s1);
	size_t len2 = strlen(s2);
	size_t len = len1 + len2;
	char *p = (char*)xmalloc(len + 1);
	memcpy(p, s1, len1);
	memcpy(p + len1, s2, len2 + 1);
	return p;
}

void get_output_filenames(const char *read_file, const char *added_suffix,
			  char **out_filename_ret, char **log_filename_ret)
{
	char *dot;
	size_t read_filename_len = strlen(read_file);
	size_t suffix_len = strlen(added_suffix);

	char read_file_copy[read_filename_len + 1];
	strcpy(read_file_copy, read_file);

	dot = strrchr(read_file_copy, '.');
	if (dot && strcmp(dot, ".gz") == 0)
		*dot = '\0';

	char *base_name = strrchr(read_file_copy, '/');
	if (base_name) {
		do {
			base_name++;
		} while (*base_name == '/');
	} else {
		base_name = read_file_copy;
	}

	size_t base_name_len = strlen(base_name);
	dot = strrchr(base_name, '.');
	if (!dot)
		dot = &base_name[base_name_len];

	char *new_read_file = (char*)xmalloc(base_name_len + 1 + suffix_len + 1);
	char *log_filename = (char*)xmalloc(base_name_len + 1 + suffix_len + sizeof(".log") + 1);

	memcpy(new_read_file, base_name, dot - base_name);
	memcpy(log_filename, base_name, dot - base_name);

	sprintf(&new_read_file[dot - base_name], ".%s%s", added_suffix, dot);
	sprintf(&log_filename[dot - base_name], ".%s.log", added_suffix);

	*out_filename_ret = new_read_file;
	*log_filename_ret = log_filename;
}
