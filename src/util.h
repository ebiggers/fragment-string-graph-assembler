#pragma once

#include <assert.h>
#include <limits.h>
#include "compiler.h"

extern void fatal_error(const char *msg, ...) __noreturn;
extern void fatal_error_with_errno(const char *msg, ...);
extern void info(const char *format, ...);

#define unimplemented() \
	fatal_error("Unimplemented at %s:%d", __FILE__, __LINE__)

#define END_LONGOPTS \
	{"help",        no_argument,       NULL, 'h'}, \
	{NULL,          0,                 NULL, 0}

#define PROCESS_OTHER_OPTS \
		default:			\
		case 'h':			\
			usage();			\
			exit(c == 'h' ? 0 : 2);

#define DEFINE_USAGE(s) \
static void usage() \
{ \
	const char *usage_str = s; \
	fputs(usage_str, stdout);  \
}

#define for_opt(c) \
	while ((c = getopt_long(argc, argv, optstring, longopts, NULL)) != -1)

#define USAGE_IF(c) if (c) {usage(); exit(2);}

extern long parse_long(const char *optstr, const char *argument,
		       long min = INT_MIN, long max = INT_MAX);

#define DIV_ROUND_UP(numerator, denominator) \
	(((numerator) + (denominator) + 1) / (denominator))
