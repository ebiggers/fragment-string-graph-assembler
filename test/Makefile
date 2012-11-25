all:out.reduced.mapped.collapsed.digraph.pdf out.reduced.mapped.collapsed.bidigraph.pdf

SAMPLE_SIZE     ?= 300
READ_LEN        ?= 100
READ_SEP        ?= 30
USE_PIRS        ?= true
MIN_OVERLAP_LEN ?= 25
GENOME          ?= random_genome.fa
BIDIGRAPH_OPS   ?= false

ifeq ($(GENOME),random_genome.fa)
genome.fa.$(SAMPLE_SIZE):
	./gen_random_genome.pl $(SAMPLE_SIZE) > $@
else
genome.fa.$(SAMPLE_SIZE): $(GENOME)
	fasta-head -c $(SAMPLE_SIZE) $+ > $@
endif

ifeq ($(USE_PIRS),true)
pirs_reads_100_180_1.fq:genome.fa.$(SAMPLE_SIZE)
	pirs simulate -x 30 --no-gc-bias --no-subst-errors --no-indels \
			--random-seed=1 --threads=1 $+
else
reads.fa:genome.fa.$(SAMPLE_SIZE) simulate_uniform_reads.pl
	./simulate_uniform_reads.pl --read-len=$(READ_LEN)	\
			--read-sep=$(READ_SEP) --allow-rc $+ > $@
endif

ifeq ($(USE_PIRS),true)
reads.bvv:pirs_reads_100_180_1.fq
	convert-reads $+ $@
else
reads.bvv:reads.fa
	convert-reads $+ $@
endif

out.overlaps:reads.bvv
	compute-overlaps $+ $@ -l $(MIN_OVERLAP_LEN)

.remove-contained-reads:reads.bvv out.overlaps
	remove-contained-reads reads.bvv reads.uncontained.bvv\
			       out.overlaps out.uncontained.overlaps \
			       out.indices_map
	touch .remove-contained-reads

reads.uncontained.bvv:.remove-contained-reads
out.uncontained.overlaps:.remove-contained-reads
out.indices_map:.remove-contained-reads

%.bidigraph.dot:%.bidigraph
	print-string-graph --dot $+ > $@

%.digraph.dot:%.digraph
	print-string-graph --dot $+ > $@

%.ps:%.dot
	dot -Tps $+ -o$@

%.pdf:%.dot
	dot -Tpdf $+ -o$@

%.png:%.dot
	dot -Tpng $+ -o$@

out.digraph:reads.uncontained.bvv out.uncontained.overlaps
	build-directed-string-graph $+ $@

out.bidigraph:reads.uncontained.bvv out.uncontained.overlaps
	build-bidirected-string-graph $+ $@

out.reduced.mapped.digraph:reads.bvv out.overlaps out.indices_map out.reduced.digraph
	map-contained-reads reads.bvv out.overlaps out.indices_map \
			    out.reduced.digraph out.reduced.mapped.digraph

%.reduced.digraph:%.digraph
	transitive-reduction $+ $@

ifeq ($(BIDIGRAPH_OPS),true)
%.reduced.bidigraph:%.bidigraph
	transitive-reduction $+ $@
endif

%.collapsed.digraph:%.digraph
	collapse-unbranched-paths $+ $@

ifeq ($(BIDIGRAPH_OPS),true)
%.collapsed.bidigraph:%.bidigraph
	collapse-unbranched-paths $+ $@
endif

ifneq ($(BIDIGRAPH_OPS),true)
%.bidigraph:%.digraph
	digraph-to-bidigraph $+ $@
endif

clean:
	rm -f reads.* out.* genome.fa.* pirs_reads* .remove-contained-reads
