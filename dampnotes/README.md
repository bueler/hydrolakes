dampnotes
==========

This is an old version of an ongoing project.  I saved it in this form at the point that I realized that the diffusive closure is actually identical to having englacial storage with a small porosity.

This project is about testing subglacial hydrology models under ice sheets, writing a paper about that, and supporting actual modeling of a tidewater glacier in Svalbard, namely Nbreen.

To build the PDF `dampnotes.pdf`, which requires LaTeX, do

    $ ln -s ~/pism/doc/ice_bib.bib   # see github.com/pism/pism
    $ make

The current model is by no means the first!  It might be the fifth or sixth in fact.  There was an early attempt to mimic the elliptic variational inequality formulation from Schoof et al. 2012; see `../phoenix/` in this repo.  (And there were even earlier models than that.)

Another model we tried is very minimal.  It works in a mode more suitable for subglacial lakes than for ice streams or surging.  It was the basis for Bueler's June 2012 talk at the Fairbanks IGS Symposium; see `bueler-igs2012/` in this repo, and the codes and data in its subdirectory `matlab/ant/`.  It has become class PISMLakesHydrology in PISM (github.com/pism/pism) and it works fine for what it is.  There is no physics for cavity or conduit evolution at all, and the water pressure is merely the overburden pressure.

The next model we tried was based on Ward van Pelt's notes of August 2012, with title "An extended hydrology model for PISM".  This model has an opening-closure equation of the type described in Hewitt (2011).  It has a damped form of the full-aquifer equation.  The damping does not obviously recover the elliptic variational inequality (in Schoof et al (2012)) in any obvious parameter limit, however.

Starting in October 2012 we changed the damping and this time we have a clear extension of the Schoof et al (2012) theory.  This model has a successful verification case and reasonable run results for Nbreen.  We understand the relationship of this model to the "lakes" model above, and to the Schoof elliptic variational inequality form, and to the Flowers & Clarke (2002) theory.  This is PISMDistributedHydrology in PISM.  The notes `dampnotes.pdf` primarily document this model.

