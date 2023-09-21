# SpecpolFlow

## NOTICE TO ALL USERS: SpecpolFlow is still under active development. Please clone the repository and update often for the latest version. The first stable release of SpecpolFlow is coming soon -- stay tuned!

SpecpolFlow is a software package that provides a completely pythonic workflow for the analysis of spectropolarimetric observations of astrophysical sources (for example, data aquired using ground-based instruments such as ESPaDOnS at CFHT, Narval at TBL, etc). It is desiged to provide a single, user-friendly pipeline from telescope to science product.

SpecpolFlow incorporates earlier pieces of software developed by Dr. Colin Folsom for the two most computationally challenging tasks: 
spectra normalization (Github: folsomcp/NormPlot) and LSD profile calculation (Github:folsomcp/LSDpy). It also provides several intermediate calculation and visualization options, including tools for developing and cleaning line masks, calculating the longitudinal magnetic field, and visualizing the LSD profile. All of these tools are fully documented in our API documentation.

We also provide and maintain a series of tutorials that can be used to teach the workflow, 
with examples of how to construct a flexible workflow from these tools for your specific needs 
(e.g., to handle automation for very large datasets using tools like pandas = Python for Data Analysis). 
These tutorials are in the form of Python notebooks, which can also be run using collaborative platforms such as Google Colab. 

The full documentation can be found here: [folsomcp.github.io/specpolFlow/](folsomcp.github.io/specpolFlow/)

Current Contributors:
* Christi Erba
* Colin Folsom
* Veronique Petit
* Shaquann Seadrow
* Patrick Stanley
* Tali Natan
* Gregg Wade
* Mary Oksala
* DJ Meleney

Past Contributors:
* Federico Villadiego Forero
* Robin Moore
* Marisol Catalan Olais
* Dax Moraes

