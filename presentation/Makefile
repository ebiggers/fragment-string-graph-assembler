presentation.pdf:presentation.tex programs.pdf
	pdflatex $<

programs.pdf:programs.dot
	dot -Tpdf -o$@ $+

clean:
	rm -f presentation.{pdf,log,aux,nav,out,snm,toc}
