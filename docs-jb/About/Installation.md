# Installation

SpecpolFlow is available as a beta release on pip:

TODO INSERT THE TEST LINK

For developer, it is also possible to install the development version from the current github repository:
```
pip install "git+https://github.com/folsomcp/LSDpy"
pip install "git+https://github.com/folsomcp/SpecpolFlow"
```

Or from any version release:
```
pip install "git+https://github.com/folsomcp/LSDpy/INSERT_RELEASE_TAG"
pip install "git+https://github.com/folsomcp/SpecpolFlow/INSERT_RELEASE_TAG"
```

## Accessing the classes and function

In interactive python (e.g. Jupyter Notebooks), importing SpecpolFlow will give access to all of the classes and function (like e.g. numpy). From these, you can build your own workflow (see the tutorials for examples).

There are two interactive tools (NormPlot and CleanMaskUI) that will not work within a notebook. These tools can be lauched from 

1. a command line 
```
Bash> python

python> import specpolflow as pol
python> my_mask = pol.read_mask('maskfile')
python> cleanMaskUI('maskfile', 'obsFile', outMaskName=None, excludeFileName='excludeRanges.dat')
```

or from a .py script containing similar instructions. 

```
> python my_mask_clean.py
```

## Using the command line scripts

the specpolFlow package also provides a few useful scripts for batch processing of files. After a pip installation, these tools (all starting with `spf-toolname`) should be available in your bash path. 

