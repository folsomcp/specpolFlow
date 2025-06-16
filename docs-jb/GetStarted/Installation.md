# Installation

SpecpolFlow is available on pip, which is the easiest way to install it:
```
pip install specpolFlow
```

There is an optional normalization tool, NormPlot, which can be installed using:
```
pip install normPlot
```
Or to install everything at once you can use:
```
pip install "specpolFlow[norm]"
```

To update SpecpolFlow and it's related packages to their most recent version use:
```
pip install --upgrade specpolFlow
pip install --upgrade LSDpy
pip install --upgrade normPlot
```

You can uninstall specpolFlow with pip:
```
pip uninstall specpolFlow
pip uninstall LSDpy
pip uninstall normPlot
```

## In-development versions

For a developer, it is also possible to install the development version from the current Github repository:
```
pip install "git+https://github.com/folsomcp/LSDpy"
pip install "git+https://github.com/folsomcp/specpolFlow"
```
This needs a relatively recent version of pip to work properly.  Version 24 or later is recommended (`pip install --upgrade pip`).


## Optional additional packages

### Packages used in tutorials

The tutorials and examples on this website largely work without installing anything extra.  However some examples illustrate using SpecpolFlow together with other packages, specifically:

* [pandas](https://pandas.pydata.org/) can be used for data analysis, and structuring or organizing data, and can be [installed with pip or conda](https://pandas.pydata.org/docs/getting_started/install.html).  There are several [optional extras](https://pandas.pydata.org/docs/getting_started/install.html#optional-dependencies) in pandas that depend on additional packages.  To get all the optional extras you can use `pip install "pandas[all]"`.  Some of our tutorials make use of fancy formatting for tables, so we recommend using at least `pip install "pandas[output-formatting]"`.  (If you installed pandas through a large distribution like Anaconda, then you likely already have these optional features.)


### Packages for developers

If you are interested in building a local copy of the documentation for SpecpolFlow, you will need to install the [Jupyter Book](https://jupyterbook.org) package.  Jupyter Book can generate webpages and PDF documents from Jupyter notebooks, and it can be [installed with pip or conda](https://jupyterbook.org/en/stable/start/overview.html).  Jupyter Book relies on [Sphinx](https://www.sphinx-doc.org/), [MyST](https://myst-parser.readthedocs.io/), and [MyST-NB](https://myst-nb.readthedocs.io/), and these should be installed automatically when you install Jupyter Book.  

## Using the command line tools

SpecpolFlow provides a several command line tools for analysis. After installing with pip, these tools (starting with `spf-toolname`) should be available in your shell's path. Some more important terminal commands are:
* `spf-makemask` generate a line mask from a line list
* `spf-cleanmask` clean problem lines from a line mask
* `spf-bz` calculate longitudinal magnetic fields
* `spf-plotlsd` plot LSD profiles
* `spf-rvfit` calculate a radial velocity from an LSD profile
* `lsdpy` calculate LSD profiles

If you are using Python in virtual environments (e.g. with [venv](https://docs.python.org/3/library/venv.html)), then you may wish to consider also installing SpecpolFlow with [pipx](https://pipx.pypa.io/stable/).  pipx will install the command line tools so that they run in their own separate virtual environment.  
```
pipx install LSDpy
pipx install specpolFlow
pipx install normPlot
```


## Accessing the classes and function

Installing SpecpolFlow provides access to the specpolFlow package in Python.  The package can be imported in other Python scripts, and in interactive Python sessions (e.g. Jupyter Notebooks):
```
import specpolFlow as pol
```
This gives access to all the classes and functions, similar to e.g. numpy.  From these, you can build your own workflow (see the [tutorials](../GetStarted/OneObservationFlow_Tutorial.ipynb) for [examples](../Tutorials/6-CalculateBz_Tutorial.ipynb)).

There are two interactive tools (NormPlot and CleanMaskUI) that will not work with remote Python sessions running on a server (e.g. Google Colab, Binder).  These usually need to be run locally, since they will open an interactive window.
