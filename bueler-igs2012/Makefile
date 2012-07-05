all: slides.pdf

slides.pdf: slides.tex
	pdflatex slides
	pdflatex slides

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

