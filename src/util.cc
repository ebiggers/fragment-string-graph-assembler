#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>

#include "util.h"

/* Prints an error message and exits the program with failure status. */
void fatal_error(const char *msg, ...)
{
	va_list va;
	va_start(va, msg);
	fflush(stdout);
	fputs("ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	putc('\n', stderr);
	va_end(va);
	exit(1);
}

void fatal_error_with_errno(const char *msg, ...)
{
	va_list va;
	int errno_save = errno;
	va_start(va, msg);
	fflush(stdout);
	fputs("ERROR: ", stderr);
	vfprintf(stderr, msg, va);
	fprintf(stderr, ": %s\n", strerror(errno_save));
	va_end(va);
	exit(1);
}

void info(const char *format, ...)
{
	va_list va;
	struct tm *tm;
	time_t t;
	static const char *month_names[] =
		{"Jan", "Feb", "Mar", "Apr", "May", "June",
		 "July", "Aug", "Sept", "Oct", "Nov", "Dec"};

	t = time(NULL);

	tm = localtime(&t);
	va_start(va, format);
	//if (info_tag)
		//fprintf(stdout, "[%s] ", info_tag);
	fprintf(stdout, "%s %d %d %02d:%02d:%02d: ",
		month_names[tm->tm_mon],
		tm->tm_mday,
		1900 + tm->tm_year,
		tm->tm_hour, tm->tm_min, tm->tm_sec);
	if (format)
		vfprintf(stdout, format, va);
	fputc('\n', stdout);
	fflush(stdout);
	va_end(va);
}

long parse_long(const char *optstr, const char *argument, long min, long max)
{
	char *tmp;
	long n = strtol(optstr, &tmp, 10);
	if (tmp == optstr || *tmp)
		fatal_error("Error parsing \"%s\": not an integer", optstr);
	if (n < min)
		fatal_error("Expected number >= %ld for argument %s", argument);
	if (n > max)
		fatal_error("Expected number <= %ld for argument %s", argument);
	return n;
}
