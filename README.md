hydrolakes
==========

This project is about constructing good subglacial hydrology models under ice sheets and glaciers, implementing and testing them in PISM, supporting actual modeling of a tidewater glacier (Nordenskioldbreen in Svalbard, i.e. "Nbreen"), and  writing a paper about that.

To build the PDF `subenhydro.pdf`, which requires LaTeX, do

    $ ln -s ~/pism/doc/ice_bib.bib   # see github.com/pism/pism
    $ make

The current model evolved from a series of more or less flawed, but mostly serious, attempts to improve the hydrology in PISM.  An old version is in `phoenix/` in this repo.  A simpler model suitable for locating subglacial lakes is documented primarily in the `bueler-igs2012/` subdirectory in this repo; in PISM it is the `-hydrology lakes` option and the PISMLakesHydrology class.  A more recent and more complete, but now out-of-date version, is in `dampnotes/` in this repo.  For the history of our attempts to build a better model see `dampnotes/README.md`.

The current model (i.e. the one described in `subenhydro.pdf`) is implemented in PISM by the `-hydrology distributed` option; this is the PISMDistributedHydrology class.  It has a successful verification case.  Running it gives reasonable results for Nbreen.

The model has an understandable relationship to:

1.  The Schoof et al (2012) distributed model, and specifically the elliptic variational inequality form of its pressure equation.
2.  The Bartholomaus et al (2011) lumped sub-/en- glacial model used for the Kennicott glacier.
3.  The Flowers & Clarke (2002) theory which assumes a function P(W), which we see is not necessary as an assumption but instead emerges as a property of steady state.

