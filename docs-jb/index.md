# Welcome to SpecpolFlow!

## About SpecpolFlow

`SpecpolFlow` is a software package that provides a completely pythonic workflow for the analysis of spectropolarimetric observations of astronomical sources (for example, data acquired using instruments such as ESPaDOnS at CFHT, Narval at TBL, etc). It is designed to provide a single, user-friendly pipeline from telescope to science product.

`SpecpolFlow`'s routines handle two computationally challenging tasks in spectropolarimetry: spectral normalization (normPlot) and least-squares deconvolution (LSD) profile calculation (LSDpy). The `SpecpolFlow` package also provides several supporting tools for developing and cleaning line masks, calculating the longitudinal magnetic field, and visualizing the LSD profile. These tools can be used through a fully documented Python API or through a command line interface.

The `SpecpolFlow` team also maintains a series of detailed tutorials with examples of how to construct a flexible workflow for your specific needs (e.g., automation for very large datasets using tools like `pandas`, the Python Data Analysis Library). These tutorials are in the form of Python notebooks, which can also be run using collaborative platforms such as Google Colab. 

## Contact us!
You can reach the SpF Development Team with questions or comments at: specpolflow@gmail.com

## Using `SpecpolFlow` in publications
To acknowledge the use of `SpecpolFlow` in publications (or talks, research results, etc.), please cite
[SpecpolFlow: a new software package for spectropolarimetry using Python](https://joss.theoj.org/papers/10.21105/joss.07891). Folsom, Erba, et al., (2025). Journal of Open Source Software, 10(111), 7891

[![DOI](https://joss.theoj.org/papers/10.21105/joss.07891/status.svg)](https://doi.org/10.21105/joss.07891)

Individual components of `SpecpolFlow` should be cited using the paper above, and can also be directly referenced as follows:
- LSDpy: The GitHub repository for LSDpy is available at https://github.com/folsomcp/LSDpy.  
- normPlot: The GitHub repository for normPlot is available at https://github.com/folsomcp/normPlot.  

A BibTeX reference for the SpecpolFlow paper is:
```
@article{SpecpolFlow2025, 
author = {{Folsom}, Colin P. and {Erba}, Christiana and {Petit}, Veronique and {Seadrow}, Shaquann and {Stanley}, Patrick and {Natan}, Tali and {Zaire}, Bonnie and {Oksala}, Mary E. and {Villadiego Forero}, Federico and {Moore}, Robin and {Catalan Olais}, Marisol},
title = {SpecpolFlow: a new software package for spectropolarimetry using Python},
journal = {Journal of Open Source Software}, 
year = {2025}, 
volume = {10}, 
number = {111}, 
pages = {7891}, 
publisher = {The Open Journal}, 
doi = {10.21105/joss.07891}, 
url = {https://doi.org/10.21105/joss.07891}
}
```