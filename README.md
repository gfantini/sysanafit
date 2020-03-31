#### SYSANAFIT

This is the sysanafit project.
It is intended as a set of tools to study decoherence in kloe data 2004-05 and related simulations.

Author:	       G. Fantini
Email:	       guido.fantini@gssi.it

Related people:	A. Di Domenico, A. De Santis
------------------------------------------------------------------------------------------------------------------
ROOT macros are tested on ROOT 6.14.06 macOS

dat/	    folder for ASCII data files (see documentation therein)
doc/ 	    folder that contains general documentation about the workflow, expected I/O of files
inc/	    libraries, general purpose functions and classes to be included in other parts of the software, theoretical model
legacy/     here some interactive macros are placed. They date back to a moment	when the core of the fit was performed by a fortran code called fitinterf
ntp/	    some data and mc simulations (ROOT and old hbook). A repository is not a storage place, please do not upload data here.
	    These files are just meant to allow running tests on the software.
src/	    source code (root macros at the moment) to perform fit



###WARNING (for posterity):
******************************************************************************************************************
At the time this code was written I was moving my first steps with C / C++ in a project bigger
than just an exercise from an undergrad programming class.
I was not even aware of how to compile a ROOT macro (this is why everything is interactive).
Also I was not aware of repositories, so I used to keep track of versioning by hand.
This is why the same file exists in different versions e.g.: file_v06.h file_v08.h
Often times such versions are not compatible with one another.

These and many other things do NOT follow the best practices and make the work of the analyst
a lot harder (and sometimes frustrating). I apologize.
I tried to make things right documenting as much as I can and putting all the work in this git repository.


### Getting started
The easiest way to get started with this version of the sw is to navigate in the src/ directory then
$ root -l
[0] .L sysanacustom_v2.0.0_test.cpp
[1] test("../ntp/Histograms_default.root");

At present all the fits will run and hopefully a canvas with the profiled likelihood will be produced.
Still there is something that causes the sw to crash.


Further data (very unorganized, practically a cut & paste of my folder where I was working up to 2017)
is available here https://mega.nz/#!topnGbaK!DlXZv4kQA_SaVQkORov04xG8sif1Yw_fJzx7Iorj5cQ