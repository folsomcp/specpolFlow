# Installation

SpecpolFlow is available on pip, which is the easy way to install it:
```
pip install specpolFlow
```

There is an optional normalization tool, NormPlot, which can be installed:
```
pip install normPlot
```
Or to install everything at once you can use:
```
pip install "specpolFlow[norm]"
```

You can also uninstall specpolFlow with pip:
```
pip uninstall specpolFlow
pip uninstall normPlot
pip uninstall LSDpy
```

## In-development versions

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

## Using the command line tools

SpecpolFlow provides a several command line tools for analysis. After installing with pip, these tools (starting with `spf-toolname`) should be available in your shell's path. Some more important terminal commands are:
* `spf-makemask` generate a line mask from a line list
* `spf-cleanmask` clean problem lines from a line mask
* `spf-bz` calculate longitudinal magnetic fields
* `spf-plotlsd` plot LSD profiles
* `spf-rvfit` calculate a radial velocity from an LSD profile
* `lsdpy` calculate LSD profiles


## Accessing the classes and function

Installing SpecpolFlow provides access to the specpolFlow package in Python.  The package can be imported in other Python scripts, and in interactive Python sessions (e.g. Jupyter Notebooks):
```
import specpolFlow as pol
```
This gives access to all the classes and functions, similar to e.g. numpy.  From these, you can build your own workflow (see the [tutorials](../GetStarted/OneObservationFlow_Tutorial.ipynb) for [examples](../Tutorials/6-CalculateBz_Tutorial.ipynb)).

There are two interactive tools (NormPlot and CleanMaskUI) that will not work with remote Python sessions running on a server (e.g. Google Colab, Binder).  These usually need to be run locally, since they will open an interactive window.
