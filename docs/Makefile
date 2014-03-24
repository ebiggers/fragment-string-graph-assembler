SHELL = bash

all:MATH469-slides.pdf MATH469-paper.pdf capstone-slides.pdf

MATH469-slides.pdf:MATH469-slides.tex programs.pdf
	pdflatex $<

capstone-slides.pdf:capstone-slides.tex programs.pdf
	pdflatex $<

programs.pdf:programs.dot
	dot -Tpdf -o$@ $+

MATH469-paper.pdf:MATH469-paper.tex refs.bib
	pdflatex $<
	bibtex MATH469-paper.aux
	pdflatex $<
	pdflatex $<


clean:
	rm -f {MATH469-paper,MATH469-slides,capstone-slides}.{pdf,log,aux,nav,out,snm,toc,bbl,blg,fff,lot,ttt,lof,vrb} programs.pdf

.PHONY: all clean
