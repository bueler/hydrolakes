all: phoenix.pdf

codes   := 
figures := verifcase1.pdf verifcase2.pdf convergecase1.pdf boundary_cartoon.png

codes   := $(addprefix matlab/, $(codes))
figures := $(addprefix figs/, $(figures))

# link ice-bib.bib needed for bibtex to work  (to file in pism-dev/doc/)

%.pdf: %.tex $(figures) ice-bib.bib
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *~

