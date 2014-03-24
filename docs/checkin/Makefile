#
# Generic Makefile for tex documents
#

SHELL := bash
TEX   := $(wildcard *.tex)
BASE  := $(basename $(TEX))
PDF   := $(addsuffix .pdf,$(BASE))
BIB   := $(wildcard *.bib)

all:$(PDF)

$(PDF): %.pdf:%.tex $(BIB)
	pdflatex $*.tex
	if grep -q '\\endnote' $*.tex; then	\
		pdflatex $*.tex;		\
	fi
	if [ "$(BIB)" != "" ]; then		\
		if grep -q 'biber' $*.tex; then \
			biber $*;		\
		else				\
			bibtex $*.aux;		\
		fi;				\
		pdflatex $*.tex;		\
		pdflatex $*.tex;		\
	fi


clean:
	for base in $(BASE); do \
		rm -f $${base}.{log,out,aux,pdf,ent,bbl,blg,snm,toc,nav,vrb,aux.blg,bcf,run.xml}; \
	done

.PHONY: all clean
