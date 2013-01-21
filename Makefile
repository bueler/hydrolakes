all: subenhydro.pdf

# link ice_bib.bib needed for bibtex to work  (to that file in pism-dev/doc/)

subenhydro.pdf: subenhydro.tex ice_bib.bib
	pdflatex subenhydro
	bibtex subenhydro
	pdflatex subenhydro
	pdflatex subenhydro

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

