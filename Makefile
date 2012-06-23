all: notes.pdf

codes   := 
figures := 

codes   := $(addprefix matlab/, $(codes))
figures := $(addprefix figs/, $(figures))

# link ice_bib.bib needed for bibtex to work  (to file in pism-dev/doc/)

%.pdf: %.tex $(figures) ice_bib.bib
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *~

