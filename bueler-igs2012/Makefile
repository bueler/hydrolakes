all: slides.pdf notes.pdf

slides.pdf: slides.tex
	pdflatex slides
	pdflatex slides

figures := diffstencil.pdf

figures := $(addprefix ../figs/, $(figures))

# link ice_bib.bib needed for bibtex to work  (to that file in pism-dev/doc/)

notes.pdf: notes.tex $(figures) ice_bib.bib
	pdflatex notes
	bibtex notes
	pdflatex notes
	pdflatex notes

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

