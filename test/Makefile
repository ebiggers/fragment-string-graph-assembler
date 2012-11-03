all:out.reduced.graph

SAMPLE_SIZE ?= 300
READ_LEN ?= 100
READ_SEP ?= 30
USE_PIRS ?= true

genome.fa.$(SAMPLE_SIZE): NC_000913.fna
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

%.dot:%.graph
	print-graph --dot $+ > $@

%.ps:%.dot
	dot -Tps $+ -o$@ -v

%.pdf:%.dot
	dot -Tpdf $+ -o$@ -v

%.png:%.dot
	dot -Tpng $+ -o$@ -v

out.graph:reads.uncontained.bvv out.uncontained.overlaps
	build-graph $+ $@

out.reduced.graph:out.graph
	transitive-reduction $+ $@

clean:
	rm -f reads.* out.* genome.fa.* pirs_reads* .remove-contained-reads
