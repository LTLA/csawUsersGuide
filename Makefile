all: csaw.pdf compile.tar.gz

compile.tar.gz: csaw.tex 
	tar -czf $@ $< plots-ug/ Bioconductor2.sty unsrturl.bst ref_ug.bib

csaw.pdf: csaw.tex
	pdflatex $<
	bibtex csaw
	pdflatex $<
	pdflatex $<

