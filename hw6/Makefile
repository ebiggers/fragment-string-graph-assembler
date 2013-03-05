#
# Generic Makefile for tex documents
#

SHELL := bash
TEX   := $(wildcard *.tex)
BASE  := $(basename $(TEX))
PDF   := $(addsuffix .pdf,$(BASE))

all:$(PDF)

$(PDF): %.pdf:%.tex
	pdflatex $+

clean:
	for base in $(BASE); do \
		rm -f $${base}.{log,out,aux,pdf}; \
	done

.PHONY: all clean
