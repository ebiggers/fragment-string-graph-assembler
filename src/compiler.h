#pragma once

#ifdef __GNUC__
#	if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
# 		define __cold __attribute__((cold))
#	else
#		define __cold
#	endif
#	define __noreturn __attribute__((noreturn))
#	define __format(type, format_str, args_start) \
			__attribute__((format(type, format_str, args_start)))
#	define __noinline __attribute__((noinline))
#else
#	define __noreturn
#	define __cold
#	define __format
#	define __noinline
#endif

/*
 * Compare-and-swap.  Equivalent to the folliwng, but executed
 * atomically:
 *
 * Q tmp = *ptr;
 * if (tmp == oval)
 * 	*ptr = nval;
 * return tmp;
 */
template <typename Q>
static inline Q cas(volatile Q *ptr, Q oval, Q nval) {
	/* gcc builtin */
	return __sync_val_compare_and_swap(ptr, oval, nval);
}

template <typename Q>
static inline bool cas_bool(volatile Q *ptr, Q oval, Q nval) {
	/* gcc builtin */
	return __sync_bool_compare_and_swap(ptr, oval, nval);
}

/*
 * Atomically writes @nval into @ptr and returns the previous value of
 * @ptr.
 */
template<typename T>
static inline T atomic_set(volatile T *ptr, T nval) {
	/* gcc builtin */
	return __sync_lock_test_and_set(ptr, nval);
}
