# Installation

SpecpolFlow is available as a beta release on pip:
```
pip install specpolFlow
```
To also install the optional spectrum normalization tool NormPlot use:
```
pip install specpolFlow[norm]
```
In some terminals, like zsh, you may need to put the command in quotes
```
pip install "specpolFlow[norm]"
```


For a developer, it is also possible to install the development version from the current Github repository:
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

In interactive python (e.g. Jupyter Notebooks), importing SpecpolFlow will give access to all of the classes and function (similar to numpy). From these, you can build your own workflow (see the tutorials for examples).

There are two interactive tools (NormPlot and CleanMaskUI) that will not work within a notebook. These tools can be launched from 

1. a command line 
```
Bash> python

python> import specpolFlow as pol
python> my_mask = pol.read_mask('maskfile.dat')
python> pol.cleanMaskUI('maskfile.dat', 'obsFile.s', outMaskName=None, excludeFileName='excludeRanges.dat')
```

2. from a .py script containing similar instructions. 

```
> python my_mask_clean.py
```

## Using the command line scripts

SpecpolFlow provides a several command line tools for analysis. After installing with pip, these tools (all starting with `spf-toolname`) should be available in your shell's path. Some more important terminal commands are:
* `spf-makemask` generate a line mask from a line list
* `spf-cleanmask` clean problem lines from a line mask
* `spf-bz` calculate longitudinal magnetic fields
* `spf-plotlsd` plot LSD profiles
* `spf-rvfit` calculate a radial velocity from an LSD profile
* `lsdpy` calculate LSD profiles
