
LDADD = libassemble.a -lboost_serialization
AM_CXXFLAGS = 

bin_PROGRAMS =				\
	bidigraph-to-digraph		\
	bidigraph-eulerian-cycle	\
	build-bidirected-string-graph	\
	build-directed-string-graph	\
	calculate-A-statistics		\
	collapse-unbranched-paths	\
	compute-overlaps		\
	convert-reads			\
	digraph-to-bidigraph		\
	extract-edge-seqs		\
	map-contained-reads		\
	min-cost-circulation		\
	print-overlaps			\
	print-string-graph		\
	remove-contained-reads		\
	transitive-reduction

noinst_LIBRARIES = libassemble.a

libassemble_a_SOURCES =			\
	AnyStringGraph.h		\
	BaseUtils.cc			\
	BaseUtils.h			\
	BaseVec.h			\
	BaseVecVec.cc			\
	BaseVecVec.h			\
	BidirectedStringGraph.cc	\
	BidirectedStringGraph.h		\
	compiler.h			\
	DirectedStringGraph.cc		\
	DirectedStringGraph.h		\
	Kmer.h				\
	Overlap.cc			\
	Overlap.h			\
	StringGraph.h			\
	util.cc				\
	util.h

bidigraph_to_digraph_SOURCES          = bidigraph-to-digraph.cc
bidigraph_eulerian_cycle_SOURCES      = bidigraph-eulerian-cycle.cc
build_bidirected_string_graph_SOURCES = build-bidirected-string-graph.cc
build_directed_string_graph_SOURCES   = build-directed-string-graph.cc
calculate_A_statistics_SOURCES        = calculate-A-statistics.cc
collapse_unbranched_paths_SOURCES     = collapse-unbranched-paths.cc
compute_overlaps_SOURCES              = compute-overlaps.cc
convert_reads_SOURCES                 = convert-reads.cc
digraph_to_bidigraph_SOURCES          = digraph-to-bidigraph.cc
extract_edge_seqs_SOURCES             = extract-edge-seqs.cc
map_contained_reads_SOURCES           = map-contained-reads.cc
min_cost_circulation_SOURCES          = min-cost-circulation.cc
print_overlaps_SOURCES                = print-overlaps.cc
print_string_graph_SOURCES            = print-string-graph.cc
remove_contained_reads_SOURCES        = remove-contained-reads.cc
transitive_reduction_SOURCES          = transitive-reduction.cc
