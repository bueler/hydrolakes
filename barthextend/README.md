hydrolakes/barthextend
==========

This correspondence will be submitted to J. Glaciology.  To build the PDF `barthextend.pdf`, which requires LaTeX, do

    $ ln -s ~/pism/doc/ice-bib.bib   # see github.com/pism/pism
    $ mkdir igs && cd igs
    $ wget http://www.igsoc.org/production/igs-v2_00-distrib.zip
    $ unzip igs-v2_00-distrib.zip
    $ cd ..
    $ ln -s igs/igs.cls
    $ ln -s igs/igs.bst
    $ ln -s igs/igsnatbib.sty
    $ make

The model described in `gmd-hydro.pdf` observes that the Bartholomaus et al (2011) lumped subglacial-englacial model used for the Kennicott glacier is a motivation (or reduction) of the main PISM `distributed` model.  This short paper documents that.

