# Generic make file for LaTeX: requires GNU make

TEXFILE	= gmd-hydro.tex

all: $(TEXFILE:.tex=.pdf)

# unzip figszip.zip to build

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

