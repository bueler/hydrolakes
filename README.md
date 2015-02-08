hydrolakes
==========

This project is about constructing good subglacial hydrology models under ice sheets and glaciers, implementing and testing them in PISM, supporting actual modeling of a tidewater glacier (Nordenskioldbreen in Svalbard, i.e. "Nbreen") and an ice sheet (Greenland), and writing a paper about that.

This paper has been submitted to Geoscientific Model Development, an EGU journal.  It is `gmd-2014-113`.

To build the PDF `gmd-hydro.pdf`, which requires LaTeX, do

    $ mkdir Copernicus
    $ cd Copernicus
    $ wget http://www.geoscientific-model-development.net/Copernicus_LaTeX_Package_v_2_7.zip
    $ unzip Copernicus_LaTeX_Package_v_2_7.zip
    $ cd ..
    $ ln -s Copernicus/copernicus.cls
    $ ln -s Copernicus/copernicus.bst
    $ make

The model described in `gmd-hydro.pdf` is implemented in PISM by the `-hydrology distributed` option; this is the `PISMDistributedHydrology` class.  It has a successful verification case.  Running it gives reasonable results for Nbreen (see van Pelt (2013)) and Greenland.  The model has an understandable relationship to several other models, but it lacks conduits; see `gmd-hydro.pdf`.

history
-------

The current model evolved from a series of more or less flawed, but mostly serious, attempts to improve the hydrology in PISM.  An old version is in `phoenix/` in this repo.  The simpler `-hydrology routing` model is observed to be suitable for locating subglacial lakes in the `bueler-igs2012/` subdirectory in this repo.  Another out-of-date version is in `dampnotes/` in this repo; for some of the history of our attempts to build a better model see `dampnotes/README.md`.  Until July 2013 the latest version was `subenhydro.tex`, and that was the model used in the last chapter of Ward's thesis.  Then the actual englacial part was removed and the name became `subhydro.tex`.  Then Andy changed it to `gmd-hydro.tex` in an attempt to get things going toward submission to _Geoscientific Model Development_.  The analysis of the pressure equation implicit in (Bartholomaus et al.~2011) was stripped out of the main paper in April 2014, and submitted by Bueler as a correspondence in _J.~Glaciol.; see subdirectory `barthextend/`.

