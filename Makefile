all: notes.pdf slides.pdf

codes   := 
figures := diffstencil.pdf enthalpy-model-crop.png hydro-mask.png siple.png ftt-mask.png pism-users-map.png

codes   := $(addprefix codes/, $(codes))
figures := $(addprefix figs/, $(figures))

# link ice_bib.bib needed for bibtex to work  (to that file in pism-dev/doc/)

notes.pdf: notes.tex $(figures) ice_bib.bib
	pdflatex notes
	bibtex notes
	pdflatex notes
	pdflatex notes

slides.pdf: slides.tex $(figures)
	pdflatex slides
	pdflatex slides

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

