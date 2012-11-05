all:out.reduced.digraph out.reduced.bidigraph

SAMPLE_SIZE ?= 300
READ_LEN ?= 100
READ_SEP ?= 30
USE_PIRS ?= true
GENOME ?= random_genome.fa

genome.fa.$(SAMPLE_SIZE): $(GENOME)
	fasta-head -c $(SAMPLE_SIZE) $+ > $@

ifeq ($(USE_PIRS),true)
pirs_reads_100_180_1.fq:genome.fa.$(SAMPLE_SIZE)
	pirs simulate -x 10 --no-gc-bias --no-subst-errors --no-indels \
			--random-seed=1 --threads=1 $+
else
reads.fa:genome.fa.$(SAMPLE_SIZE) simulate_uniform_reads.pl
	./simulate_uniform_reads.pl --read-len=$(READ_LEN)	\
			--read-sep=$(READ_SEP) $+ > $@
endif

ifeq ($(USE_PIRS),true)
reads.bvv:pirs_reads_100_180_1.fq
	convert-reads $+ $@
else
reads.bvv:reads.fa
	convert-reads $+ $@
endif

out.overlaps:reads.bvv
	compute-overlaps $+ $@

.remove-contained-reads:reads.bvv out.overlaps
	remove-contained-reads reads.bvv reads.uncontained.bvv\
			       out.overlaps out.uncontained.overlaps
	touch .remove-contained-reads

reads.uncontained.bvv:.remove-contained-reads

out.uncontained.overlaps:.remove-contained-reads

%.bidigraph.dot:%.bidigraph
	print-string-graph --dot $+ > $@

%.digraph.dot:%.digraph
	print-string-graph --dot $+ > $@

%.ps:%.dot
	dot -Tps $+ -o$@ -v

%.pdf:%.dot
	dot -Tpdf $+ -o$@ -v

%.png:%.dot
	dot -Tpng $+ -o$@ -v

out.digraph:reads.uncontained.bvv out.uncontained.overlaps
	build-directed-string-graph $+ $@

out.bidigraph:reads.uncontained.bvv out.uncontained.overlaps
	build-bidirected-string-graph $+ $@

%.reduced.digraph:%.digraph
	transitive-reduction $+ $@

%.reduced.bidigraph:%.bidigraph
	transitive-reduction $+ $@

clean:
	rm -f reads.* out.* genome.fa.* pirs_reads* .remove-contained-reads
