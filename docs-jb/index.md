# Welcome to SpecpolFlow

:::{warning}
The SpecpolFlow documentation is currently undergoing a major restructuration. 
Stay tuned!
:::

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

:::{Warning}
NOTICE TO ALL USERS: SpecpolFlow is still under active development. Please clone the repository and update often for the latest version. The first stable release of SpecpolFlow is coming soon -- stay tuned!
:::