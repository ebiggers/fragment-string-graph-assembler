SHELL = bash

all:presentation.pdf paper.pdf

presentation.pdf:presentation.tex programs.pdf
	pdflatex $<

programs.pdf:programs.dot
	dot -Tpdf -o$@ $+

paper.pdf:paper.tex refs.bib
	pdflatex $<
	bibtex paper.aux
	pdflatex $<
	pdflatex $<


clean:
	rm -f {paper,presentation}.{pdf,log,aux,nav,out,snm,toc,bbl,blg,fff,lot,ttt,lof,vrb}

.PHONY: all clean
