hydrolakes
==========

This project is about constructing good subglacial hydrology models under ice sheets and glaciers, implementing and testing them in PISM, supporting actual modeling of a tidewater glacier (Nordenskioldbreen in Svalbard, i.e. "Nbreen"), and  writing a paper about that.

This will be submitted to Geoscientific Model Development, an EGU journal.  To build the PDF `gmd-hydro.pdf`, which requires LaTeX, do

    $ ln -s ~/pism/doc/ice-bib.bib   # see github.com/pism/pism
    $ mkdir Copernicus
    $ cd Copernicus
    $ wget http://www.geoscientific-model-development.net/Copernicus_LaTeX_Package_v_2_7.zip
    $ unzip Copernicus_LaTeX_Package_v_2_7.zip
    $ cd ..
    $ ln -s Copernicus/copernicus.cls
    $ ln -s Copernicus/copernicus.bst
    $ make

The model described in `gmd-hydro.pdf` is implemented in PISM by the `-hydrology distributed` option; this is the `PISMDistributedHydrology` class.  It has a successful verification case.  Running it gives reasonable results for Nbreen.

The model has an understandable relationship to several other models:

1.  The model in which subglacial water pressure equals overburden pressure, i.e. Shreve (1972).  This model is widely used for locating subglacial lakes, e.g. Siegert (2009).  This model is implemented in PISM by the `-hydrology routing` option and the `PISMRoutingHydrology` class.
2.  The Schoof et al (2012) distributed model, and specifically the elliptic variational inequality form of its pressure equation.
3.  The Bueler & Brown (2009) till-only model where there is a reasonably-well-understood relation between a highly-simplified 'pore water' pressure and a till yield stress.  This model is implemented in PISM by the `-hydrology null` option and the `PISMNullTransportHydrology` class.
4.  The Bartholomaus et al (2011) lumped sub-/en- glacial model used for the Kennicott glacier.
5.  The Flowers & Clarke (2002) theory which assumes a function P(W), which we see is not necessary as an assumption but instead emerges as a property of steady state.

But the model lacks conduits.

The current model evolved from a series of more or less flawed, but mostly serious, attempts to improve the hydrology in PISM.  An old version is in `phoenix/` in this repo.  The simpler `-hydrology routing` model is observed to be suitable for locating subglacial lakes in the `bueler-igs2012/` subdirectory in this repo.  Another out-of-date version is in `dampnotes/` in this repo; for some of the history of our attempts to build a better model see `dampnotes/README.md`.  Until July 2013 the latest version was `subenhydro.tex`, and that was the model used in the last chapter of Ward's thesis, but then the englacial part was removed, the name because `subhydro.tex`, and then Andy changed it to `gmd-hydro.tex` in an attempt to get things going toward submission to _Geoscientific Model Development_.  The analysis of the pressure equation implicit in (Bartholomaus et al.~2011) was stripped out of the main paper in April 2014, for submission by Bueler as a correspondence in _J.~Glaciol._


