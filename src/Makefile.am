
LDADD = libassemble.a -lboost_serialization
AM_CXXFLAGS = -std=c++11

bin_PROGRAMS =			\
	convert-reads		\
	compute-overlaps	\
	print-overlaps		\
	remove-contained-reads	\
	build-graph		\
	print-graph		\
	transitive-reduction

noinst_LIBRARIES = libassemble.a

libassemble_a_SOURCES =  \
	Base.h		\
	BaseUtils.cc	\
	BaseUtils.h	\
	BaseVec.h	\
	BaseVecVec.cc	\
	BaseVecVec.h	\
	DirectedStringGraph.h \
	DirectedStringGraph.cc \
	BidirectedStringGraph.h \
	BidirectedStringGraph.cc \
	Overlap.cc	\
	compiler.h	\
	util.cc		\
	util.h

convert_reads_SOURCES          = convert-reads.cc
compute_overlaps_SOURCES       = compute-overlaps.cc
print_overlaps_SOURCES         = print-overlaps.cc
remove_contained_reads_SOURCES = remove-contained-reads.cc
build_graph_SOURCES            = build-graph.cc
print_graph_SOURCES            = print-graph.cc
transitive_reduction_SOURCES   = transitive-reduction.cc
