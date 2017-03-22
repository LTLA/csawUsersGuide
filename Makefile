all: csaw.pdf output.pdf

csaw.pdf: csaw.tex
	pdflatex $<
	bibtex csaw
	pdflatex $<
	pdflatex $<

output.pdf: csaw.pdf
	Rscript transferUG.R $< $@
