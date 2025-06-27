# Welcome to SpecpolFlow!

## NOTICE TO ALL USERS: 
SpecpolFlow's beta release is live at:
```
pip install specpolFlow
```
The code is still in active development, while the current modules are well tested, some details may change.

Installing SpecpolFlow will automatically install [LSDpy](https://github.com/folsomcp/LSDpy/).

The spectrum normalization tool [NormPlot](https://github.com/folsomcp/normPlot/) is an optional extra.  It can be installed at the same time as SpecpolFlow using:
```
pip install specpolFlow[norm]
```
or NormPlot can be installed separately using:
```
pip install normPlot
```


## About SpecpolFlow

SpecpolFlow is a software package that provides a completely pythonic workflow for the analysis of spectropolarimetric observations of stellar sources (for example, data acquired from ESPaDOnS at CFHT, Narval at TBL, etc). It is designed to provide a single, user-friendly pipeline from telescope to science product.

SpecpolFlow incorporates tools for spectra normalization (Github: [folsomcp/normPlot](https://github.com/folsomcp/normPlot)) and LSD profile calculation (Github: [folsomcp/LSDpy](https://github.com/folsomcp/LSDpy)). It also provides several supporting tools, including tools for developing and cleaning line masks, calculating the longitudinal magnetic field, and visualizing LSD profiles. These tools can be used through a fully documented Python API, or through a command line interface.

We also maintain a series of detailed tutorials, with examples of how to construct a flexible workflow from SpecpolFlow's tools for your specific needs (e.g., automation for very large datasets using tools like pandas = Python for Data Analysis). 
These tutorials are in the form of Jupyter Notebooks, which can also be run using collaborative platforms such as Google Colab. 

The full documentation can be found here: [folsomcp.github.io/specpolFlow/](https://folsomcp.github.io/specpolFlow/)

## Contact us!
You can reach the SpF Development Team at: specpolflow@gmail.com

## SpF Development Team:
* Christi Erba (co-PI)
* Colin Folsom (co-PI)
* Veronique Petit
* Shaquann Seadrow
* Patrick Stanley
* Tali Natan
* Bonnie Zaire

Current Contributors:
* Gregg Wade
* Mary Oksala

Past Contributors:
* Federico Villadiego Forero
* DJ Meleney
* Robin Moore
* Dax Moraes
* Marisol Catalan Olais

## Logo:
The SpecpolFlow logo was created by the talented Tali Natan! Please contact the SpF team if you would like to use this graphic in a publication. 

