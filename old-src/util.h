/*
 * util.h
 *
 * Miscellaneous error, warning, file, memory, and threading functions.
 */

#ifndef _UTIL_H
#define _UTIL_H

#include <stdio.h>
#include <stddef.h>
#include <sys/types.h>
#include "compiler.h"
#include <ctype.h>

enum OutputType {
	TEXT,
	GZIP,
	DEFAULT,
};

extern unsigned get_default_num_threads() __cold;
extern void fatal_error(const char *msg, ...) __noreturn __cold __format(printf, 1, 2);
extern void fatal_error_with_errno(const char *msg, ...) __noreturn __cold __format(printf, 1, 2);
extern void warning(const char *msg, ...) __cold __format(printf, 1, 2);
extern void info(const char *format, ...) __cold __format(printf, 1, 2);
extern const char *info_tag;
extern char *xstrdup(const char *s);
extern char *xstrdup2(const char *s1, const char *s2);
extern void *xmalloc(size_t size);
extern void *xzalloc(size_t size);
extern void *xrealloc(void *ptr, size_t size);
extern void *xmemalign(size_t alignment, size_t size);

extern void get_output_filenames(const char *read_file,
				 const char *added_suffix,
				 char **out_filename_ret,
				 char **log_filename_ret);

inline void __noreturn __cold __unimplemented(const char *filename, int line)
{
	fatal_error("%s:%d: Unimplemented code", filename, line);
}

/* Remove all whitespace from the end of the line/string.  Return the length of
 * the trimmed string. */
static inline size_t trim(char *s, size_t len)
{
	while (len != 0 && isspace(s[len - 1]))
		s[--len] = '\0';
	return len;
}

#define unimplemented() __unimplemented(__FILE__, __LINE__)

#endif
