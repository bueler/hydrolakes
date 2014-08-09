# Generic make file for LaTeX: requires GNU make

TEXFILE	= barthextend.tex

all: $(TEXFILE:.tex=.pdf)

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

