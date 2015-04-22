hydrolakes
==========

This project is about constructing good subglacial hydrology models for ice sheets and glaciers, implementing and testing them in PISM, supporting actual modeling of a tidewater glacier (Nordenskioldbreen in Svalbard) and an ice sheet (Greenland), and writing a paper about that.

This paper has been submitted and reviewed in Geoscientific Model Development Discussions, an EGU journal.  It is [`gmd-2014-113`](http://www.geosci-model-dev-discuss.net/7/4705/2014/).

The model described in the paper, which has both till and a linked cavity network, is run in PISM using the `-hydrology distributed` option; this is the `PISMDistributedHydrology` class in the [source code](https://github.com/pism/pism/tree/dev/src/base/hydrology).  A simplification without linked cavities, which instead uses overburden pressure as the modeled water pressure, uses option `-hydrology routing`; it is class `PISMRoutingHydrology`.

The model has a successful verification case.  Running it gives reasonable results for Nordenskioldbreen (see [van Pelt (2013)](http://wardvanpelt.com/thesis_WvanPelt.pdf)) and Greenland (the [paper](http://www.geosci-model-dev-discuss.net/7/4705/2014/)).  The model has an understandable relationship to several other models, but it lacks conduits.

Using the model
---------------

The model is part of [PISM](http://www.pism-docs.org).  It appeared in the `stable0.6` release (February 2014) and will be supported, and possibly improved, in future versions.  The source code is in directory `src/base/hydrology/` in the PISM download.

Build the PDF
-------------

To build `gmd-hydro.pdf`, which requires LaTeX and the Copernicus LaTeX class, do

    $ mkdir Copernicus
    $ cd Copernicus
    $ wget http://www.geoscientific-model-development.net/Copernicus_LaTeX_Package_v_2_7.zip
    $ unzip Copernicus_LaTeX_Package_v_2_7.zip
    $ cd ..
    $ ln -s Copernicus/copernicus.cls
    $ ln -s Copernicus/copernicus.bst
    $ make

History
-------

The current model evolved from a series of more or less flawed, but mostly serious, attempts to improve the hydrology in PISM.  Note that work started in summer 2012.

An old version is in `phoenix/` in this repo.  The simpler `-hydrology routing` model is observed to be suitable for locating subglacial lakes in the `bueler-igs2012/` subdirectory in this repo.  Another out-of-date version is in `dampnotes/` in this repo; for some of the history of our attempts to build a better model see `dampnotes/README.md`.  Until July 2013 the latest version was `subenhydro.tex`, and that was the model used in the last chapter of Ward's thesis.  Then the actual englacial part was removed and the name became `subhydro.tex`.  Then Andy changed it to `gmd-hydro.tex` in an attempt to get things going toward submission to _Geoscientific Model Development_.

An analysis of the pressure equation implicit in (Bartholomaus et al.~2011) was stripped out of the main paper in April 2014, and published in the [_Journal of Glaciology_](http://www.ingentaconnect.com/content/igsoc/jog/2014/00000060/00000222/art00018).  See subdirectory `barthextend/`.

