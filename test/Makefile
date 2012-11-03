all:out.graph

SAMPLE_SIZE = 5000

genome.fa.$(SAMPLE_SIZE): NC_000913.fna
	fasta-head -c $(SAMPLE_SIZE) $+ > $@

pirs_reads_100_180_1.fq:genome.fa.$(SAMPLE_SIZE)
	pirs simulate -x 20 --no-gc-bias --no-subst-errors --no-indels --no-log-files \
			--random-seed=1 --threads=1 $+

reads.bvv:pirs_reads_100_180_1.fq
	convert-reads $+ $@

out.overlaps:reads.bvv
	compute-overlaps $+ $@

out.dot:out.graph
	print-graph --dot $+ > $@

#out.ps:out.dot
	#dot -Tps $+ -o$@

#out.pdf:out.dot
	#dot -Tpdf $+ -o$@

out.graph:reads.bvv out.overlaps
	build-graph $+ $@

clean:
	rm -f reads.* out.* genome.fa.* pirs_reads*
