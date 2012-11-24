#pragma once

#include <assert.h>
#include <limits.h>
#include "compiler.h"
#include <stdio.h>


#if __cplusplus >= 201103L
#define foreach(member, list) \
	for (member : list)

#else
#include <boost/foreach.hpp>
#define foreach(member, list) \
		BOOST_FOREACH(member, list)
#endif


extern void fatal_error(const char *msg, ...) __noreturn;
extern void fatal_error_with_errno(const char *msg, ...);
extern void info(const char *format, ...);

#define unimplemented() \
	fatal_error("Unimplemented at %s:%d", __FILE__, __LINE__)

#define unreachable() \
	fatal_error("unreachable() at %s:%d", __FILE__, __LINE__)

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

#define USAGE_IF(cond) if (cond) {usage(); exit(2);}

extern long long parse_long(const char *optstr, const char *argument,
			    long long min = INT_MIN, long long max = INT_MAX);

#define assert2 assert


#define DIV_ROUND_UP(numerator, denominator) \
	(((numerator) + (denominator) - 1) / (denominator))

#define MODULO_NONZERO(numerator, denominator) \
	(((numerator) % (denominator)) ? ((numerator) % (denominator)) : (denominator))

#define DIV_NONZERO(numerator, denominator) \
	(((denominator) == 0) ? 0 : ((numerator) / (denominator)))

#define DOUBLE_DIV_NONZERO(numerator, denominator) \
	DIV_NONZERO(double(numerator), double(denominator))

#define TO_PERCENT(numerator, denominator) \
	(100.0 * DOUBLE_DIV_NONZERO(numerator, denominator))
