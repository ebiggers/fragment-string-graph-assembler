#pragma once

#include "StringGraph.h"
#include "DirectedStringGraph.h"
#include "BidirectedStringGraph.h"

class Overlap;

//
// Operate on either a BidirectedStringGraph or a DirectedStringGraph.
// This probably should be a class with virtual methods...
//
class AnyStringGraph {
private:
	void * impl;
	bool is_bidigraph;
public:
	AnyStringGraph(const char * filename)
	{
		impl = NULL;
		try {
			impl = (void*)new BidirectedStringGraph(filename);
			is_bidigraph = true;
		} catch(std::exception e) {
			impl = (void*)new DirectedStringGraph(filename);
			is_bidigraph = false;
		}
		if (!impl) {
			fatal_error("Error reading string graph from \"%s\"",
				    filename);
		}
	}
	~AnyStringGraph()
	{
		if (impl) {
			if (is_bidigraph)
				delete static_cast<BidirectedStringGraph*>(impl);
			else
				delete static_cast<DirectedStringGraph*>(impl);
		}
	}
#define DISPATCH0(f) \
	void f() \
	{ \
		if (is_bidigraph) \
			static_cast<BidirectedStringGraph*>(impl)->f(); \
		else	\
			static_cast<DirectedStringGraph*>(impl)->f(); \
	}
#define DISPATCH1(f, type1) \
	void f(type1 arg1) \
	{ \
		if (is_bidigraph) \
			static_cast<BidirectedStringGraph*>(impl)->f(arg1); \
		else	\
			static_cast<DirectedStringGraph*>(impl)->f(arg1); \
	}
#define DISPATCH2(f, type1, type2) \
	void f(type1 arg1, type2 arg2) \
	{ \
		if (is_bidigraph) \
			static_cast<BidirectedStringGraph*>(impl)->f(arg1, arg2); \
		else	\
			static_cast<DirectedStringGraph*>(impl)->f(arg1, arg2); \
	}
#define DISPATCH3(f, type1, type2, type3) \
	void f(type1 arg1, type2 arg2, type3 arg3) \
	{ \
		if (is_bidigraph) \
			static_cast<BidirectedStringGraph*>(impl)->f(arg1, arg2, arg3); \
		else	\
			static_cast<DirectedStringGraph*>(impl)->f(arg1, arg2, arg3); \
	}

	DISPATCH0(transitive_reduction);
	DISPATCH0(collapse_unbranched_paths);
	DISPATCH1(print_stats, std::ostream &);
	DISPATCH1(write, const char *);
	DISPATCH2(print, std::ostream &, bool);
	DISPATCH2(print_dot, std::ostream &, bool);
	DISPATCH3(map_contained_read, size_t, const Overlap &, size_t);
};
