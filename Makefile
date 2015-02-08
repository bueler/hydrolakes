# make file for LaTeX: requires GNU make

figures = Ntilfunctions.pdf exact-P-plot.pdf exact-W-plot-onu.pdf diffstencil.pdf refineWPpism.pdf g2km-init-bmelt.png g2km-init-velbase-mag.png routing-decoupled-tillwat.png routing-decoupled-bwat.png detail-routing-decoupled-bwat.png distributed-decoupled-bwat.png distributed-decoupled-bwprel.png bin1-g2km.png bin10-g2km.png bin30-g2km.png bin100-g2km.png bin300-g2km.png psteady-Po.pdf psteady-vb.pdf

all: gmd-hydro.pdf

# unzip figszip.zip to build
figszip.zip:
	mkdir figszip/
	(cd figs/ && cp ${figures} ../figszip/)
	zip figszip.zip figszip/*
	rm -rf figszip/

%.pdf: %.tex figszip.zip
	unzip figszip.zip
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~
	@rm -rf figszip/

