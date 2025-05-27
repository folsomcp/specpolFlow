# Welcome to SpecpolFlow!

## About SpecpolFlow

`SpecpolFlow` is a software package that provides a completely pythonic workflow for the analysis of spectropolarimetric observations of astronomical sources (for example, data aquired using instruments such as ESPaDOnS at CFHT, Narval at TBL, etc). It is designed to provide a single, user-friendly pipeline from telescope to science product.

`SpecpolFlow`'s routines handle two computationally challenging tasks in spectropolarimetry: spectral normalization (normPlot) and least-squares deconvolution (LSD) profile calculation (LSDpy). The `SpecpolFlow` package also provides several supporting tools for developing and cleaning line masks, calculating the longitudinal magnetic field, and visualizing the LSD profile. These tools can be used through a fully documented Python API or through a command line interface.

The `SpecpolFlow` team also maintains a series of detailed tutorials with examples of how to construct a flexible workflow for your specific needs (e.g., automation for very large datasets using tools like `pandas`, the Python Data Analysis Library). These tutorials are in the form of Python notebooks, which can also be run using collaborative platforms such as Google Colab. 

## Contact us!
You can reach the SpF Development Team with questions or comments at: specpolflow@gmail.com

## Using `SpecpolFlow` in publications
To acknowledge the use of `SpecpolFLow` in publications (or talks, research results, etc.), please cite  
&nbsp;&nbsp;&nbsp;[SpecpolFlow: a new software package for spectropolarimetry using Python.](https://arxiv.org/abs/2505.18476) Folsom, Erba, et al. (2025) Journal of Open Source Software (in review).

Individual components of `SpecpolFLow` can be referenced as follows:  
&nbsp;&nbsp;&nbsp;LSDpy: The GitHub repository for LSDpy is available at https://github.com/folsomcp/LSDpy.  
&nbsp;&nbsp;&nbsp;normPlot: The GitHub repository for normPlot is available at https://github.com/folsomcp/normPlot.  

