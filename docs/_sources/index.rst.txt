.. specpol-flow documentation master file, created by
   sphinx-quickstart on Mon Jun  7 13:58:56 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SpecPolFlow
========================================

NOTE: This documentation has been build for the RefactoringProject development branch. 

The goal of SpecpolFlow is to provide a completely pythonic workflow for the analysis of 
spectropolarimetric observations (for example, ESPaDOnS at CFHT, Narval at TBL, etc). 

The SpecpolFlow package is built around two other packages developed by Colin Folsom for 
the two most computationally challenging tasks: spectra normalization (Github: folsomcp/NormPlot) 
and LSD profile calculations (Github:folsomcp/LSDpy). 

SpecPolFlow provides intermediate calculation tools (for example line mask calculations and 
longitudinal field calculations), as well as object classes and manipulation tools that are useful 
to build a spectropolarimetric workflow (for example a LSD object class, with built-in read, write, 
slice, plot, etc methods). 

These tools are fully documented (see API documentation)

We also provide tutorials that can be used to teach the workflow, 
and examples on how to construct a workflow from these tools to e.g. handle automation 
for very large datasets (using tools like the Panda package that is becoming very popular in data science). 

These tutorials are in the form of python notebooks (that can also run on Google Colab). 


.. toctree::
   Installation
   Credits
   :maxdepth: 1
   :caption: Contents:



.. toctree::
   tutorials/OneObservationFlow_Tutorial
   tutorials/LoopFlow_Tutorial
   tutorials/CalculateBz_Tutorial
   tutorials/ConvertToSFiles_Tutorial
   tutorials/ExcludeMaskRegionClass_Tutorial
   tutorials/LSDClass_Tutorial
   tutorials/MaskClass_Tutorial
   tutorials/ColabSetup_Tutorial
   tutorials/LSDsonification_Tutorial
   :maxdepth: 1
   :caption: Tutorials:



.. toctree::
   API
   :maxdepth: 4
   :caption: API documentation:





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
