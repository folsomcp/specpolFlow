.. SpecpolFlow documentation master file, created by
   sphinx-quickstart on Thu Sep  7 12:37:15 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SpecpolFlow!
=======================================

SpecpolFlow is a software package that provides a completely pythonic workflow for the analysis of spectropolarimetric observations of astrophysical sources
(for example, data aquired using ground-based instruments such as ESPaDOnS at CFHT, Narval at TBL, etc). It is desiged to provide a single, user-friendly pipeline 
from telescope to science product.

SpecpolFlow incorporates earlier pieces of software developed by Dr. Colin Folsom for the two most computationally challenging tasks: 
spectra normalization (Github: folsomcp/NormPlot) and LSD profile calculation (Github:folsomcp/LSDpy). It also provides several intermediate 
calculation and visualization options, including tools for developing and cleaning line masks, calculating the longitudinal magnetic field, and
visualizing the LSD profile. All of these tools are fully documented in our API documentation.

We also provide and maintain a series of tutorials that can be used to teach the workflow, 
with examples of how to construct a flexible workflow from these tools for your specific needs 
(e.g., to handle automation for very large datasets using tools like pandas = Python for Data Analysis). 
These tutorials are in the form of Python notebooks, which can also be run using collaborative platforms such as Google Colab. 

.. toctree::
   Installation
   Credits
   :maxdepth: 1
   :caption: Contents

.. toctree::
   tutorials/CalculateBz
   tutorials/ColabSetup
   tutorials/LSDsonification
   tutorials/LSDtutorial
   :maxdepth: 1
   :caption: Tutorials

.. toctree::
   API
   :maxdepth: 2
   :caption: API documentation

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
