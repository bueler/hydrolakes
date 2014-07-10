# Generic make file for LaTeX: requires GNU make

TEXFILE	= gmd-hydro.tex

all: $(TEXFILE:.tex=.pdf)

# unzip figszip.zip to build

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

# link ice-bib.bib needed for bibtex to work  (to that file in pism-dev/doc/)
#%.pdf: %.tex ice-bib.bib
#	pdflatex $<
#	bibtex $(<:.tex=.aux)
#	pdflatex $<
#	pdflatex $<

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

