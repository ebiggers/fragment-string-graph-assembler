#pragma once

#include "StringGraph.h"
#include "DirectedStringGraph.h"
#include "BidirectedStringGraph.h"

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

	DISPATCH0(transitive_reduction);
	DISPATCH0(collapse_unbranched_paths);
	DISPATCH1(print_stats, std::ostream &);
	DISPATCH1(print, std::ostream &);
	DISPATCH1(print_dot, std::ostream &);
	DISPATCH1(write, const char *);
};
