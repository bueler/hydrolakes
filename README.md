hydrolakes
==========

This small project is about testing subglacial hydrology models under ice sheets.

The first model we try works in a mode more suitable for subglacial lakes than for ice streams or surging.  This is described in `notes.pdf` and in the codes `matlab/conserve.m` and `matlab/antlakes.m`, based on Anne Le Brocq's Antarctica.  It is from June 2012 and was the basis for Bueler's talk at IGS 2012; see slides in `bueler-igs2012/`.

The second model is based on Ward van Pelt's notes of August 2012, with title "An extended hydrology model for PISM", but also it relates to the earlier joint attempt in directory `phoenix/`.  This second model, like the earlier attempt, has an opening-closure equation but this time it is explicit through a damped form of the full-aquifer equation.  This new model is described in PDF `dampnotes.pdf`.

To build the two PDFs, which require LaTeX, do

    $ ln -s ~/pism/doc/ice_bib.bib
    $ make

