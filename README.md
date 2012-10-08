hydrolakes
==========

This project is about testing subglacial hydrology models under ice sheets.

To build the PDF `dampnotes.pdf`, which requires LaTeX, do

    $ ln -s ~/pism/doc/ice_bib.bib   # see github.com/pism/pism
    $ make

The current model is by no means the first!  It might be the fifth or sixth in fact.

There was an early attempt to mimic the elliptic variational inequality formulation from Schoof et al. 2012; see `phoenix/`.  (But there were even earlier models than that.)

Another model we tried is very minimal.  It works in a mode more suitable for subglacial lakes than for ice streams or surging.  It was the basis for Bueler's June 2012 talk at the Fairbanks IGS Symposium; see `bueler-igs2012/` and the codes and data in `matlab/ant/`.

The next model we tried was based on Ward van Pelt's notes of August 2012, with title "An extended hydrology model for PISM".  This model has an opening-closure equation.  It has a damped form of the full-aquifer equation.  The damping does not obviously recover the elliptic variational inequality in any obvious parameter limit, however.  This new model, or at least a close relative, is described as an alternative `dampnotes.pdf`.
