all: dampnotes.pdf

figures := diffstencil.pdf

figures := $(addprefix figs/, $(figures))

# link ice_bib.bib needed for bibtex to work  (to that file in pism-dev/doc/)

dampnotes.pdf: dampnotes.tex $(figures) ice_bib.bib
	pdflatex dampnotes
	bibtex dampnotes
	pdflatex dampnotes
	pdflatex dampnotes

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

