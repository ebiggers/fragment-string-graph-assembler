#pragma once

#include <assert.h>
#include "compiler.h"

extern void fatal_error(const char *msg, ...) __noreturn;
extern void fatal_error_with_errno(const char *msg, ...);
extern void info(const char *format, ...);

#define unimplemented() \
	fatal_error("Unimplemented at %s:%d", __FILE__, __LINE__)
