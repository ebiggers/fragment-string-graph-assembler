all:out.graph

reads.bvv:pirs_reads_100_180_1.fq
	convert-reads $+ $@

out.overlaps:reads.bvv
	compute-overlaps $+ $@

out.graph:reads.bvv out.overlaps
	build-graph $+ $@

clean:
	rm -f reads.* out.*
