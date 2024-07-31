# How to Normalize Echelle Spectra with normPlot

:::{Warning}
Still under construction. Some commands and other details may be out of date.
:::

## Introduction

This tutorial introduces the procedure of normalizing echelle spectra and the use of normPlot. The interactive normalization GUI is a part of the SpecpolFlow pipeline as a step that proceeds the polarimetric analysis. In this introduction to spectra normalization, we will be using ESPaDOnS observations of the star $\xi^1$ CMa. These are the same data that are used in the LSD and Mean Longitudinal Field tutorials, so skipping this tutorial will not prevent you from advancing further in the pipeline (if you'd like, you're welcome to use normalized spectra you create while completing this tutorial and see how the magnetic field measurements compare). 

:::{note}
The current version of this tutorial focuses on executing normPlot from the command line.
You can also run normPlot from within a Jupyter Notebook or other Python script with
```
import normPlot
normPlot.normplot('[observation_file]')
```
:::


## Converting from the .fits format to the .s text file

::::{margin}
:::{note}
ESPaDOns `p.fits` files contain a main data table with 24 columns. Here we are mostly interested in columns 19-24, which contain the unnormalized spectrum (without a telluric velocity correction), and have columns for wavelength, intensity, polarization, null 1, null 2, and errorbars.
:::
:::: 

Data files downloaded from the CADC archive are distributed in a '.fits' format. To use normPlot, these will need to be converted to the standard '.s' text file. Note that file names ending in `i.fits` contain spectroscopic data (usually 4 spectra for each polarimetric observation), and file names ending in `p.fits` contain the spectopolarimetric data. We want to convert the `p.fits` files. 

For this tutorial, we will use the ESPaDOns observation `2378196p.fits`, so the example is specific to ESPaDOnS observations from the standard Upena/LibreESPRIT pipeline at CFHT. For observations from other instruments or other data reduction pipelines, you will need a different script to convert FITS files to the `.s` format. 

::::{margin}
:::{admonition} Converting File Formats
:class: attention
See also the [tutorial on converting to .s files](../Tutorials/1-ConvertToSFiles_Tutorial.ipynb) for more information. Tools have been provided for a few commonly used spectropolarimeters (documented in the [API](https://folsomcp.github.io/specpolFlow/API/Converters_API.html)), but in some cases, you may need to write your own conversion script.
:::
::::

In the command line, write the command: 

 ``` 
 spf-fitstos-espadons inputfile-p.fits
 ```

replacing the `inputfile-p.fits` with the name of our ESPaDOnS file `2378196p.fits`. The script will return `2378196pn.s` (normalized data) and `2378196pu.s` (unnormalized data) in the same directory as the '.fits' file.

::::{margin}
:::{admonition} For more details...
:class: seealso
More details about the difference between `n.s` and `u.s` files can be found in [this tutorial](../Tutorials/1-ConvertToSFiles_Tutorial.ipynb) and in the [API](https://folsomcp.github.io/specpolFlow/API/Converters_API.html).
:::
::::

## Normalization GUI 

**Our goal is to fit a low order polynomial through "good" continuum points in every spectral order. The continuum can then be normalized by dividing the observation by the fit. This creates a ‘common line’ (a flat horizontal line at y = 1) which can be used to consistently measure the properties of the spectral lines.** 

Open the normalization GUI from the command line with
```
normplot <observation_file.s>
```

This generates a separate window with the spectrum separated into orders, which are indicated by different colors. The polynomial curves are fits to the continuum for each order, colored differently for each order. The black points are the points the used to generate the polynomial fit.

```{image} ../normplot_images/gui_open.png
:alt: gui_open
:class: bg-primary mb-1
:width: 600px
:align: center
```

The bottom left corner of the window contains a toolbar with commands from matplotlib (panning, zooming, returning to default, saving a figure). Above that are a set of buttons with some custom commands for automatically scaling both axes, auto-scaling just the y-axis, zooming, and panning. There are keyboard shortcuts for these commands (arrow keys for panning, 'a' for auto-scaling, 'A' for auto-scaling the y-axis, 'z' for zooming on a selected region, 'i' and 'o' for zooming in and out).

### The fit points
We want to select a relatively small number of "good" continuum points which can be used to generate a polynomial fit describing the continuum of the observation. These "good" points, which appear as black dots on the plot, are selected within an order and are used to fit a polynomial function for the continuum for that order. Ideally, the resulting normalized continuum will be equal to 1.

To find "good" points, each order is broken down into a set of consecutive search bins with respect to wavelength. The bins are equal in width (in velocity). In each bin, the highest flux point is selected as being most likely in the continuum (i.e. not in an absorption line). However, to deal with noise in the observation, a running average is applied to the spectrum before selecting the highest point in the search bin. This averaging helps ensure that the selected highest point in a bin is not highest because of a peak in the noise, but instead highest because of real features in the spectrum. 

## Normalizing the Individual Orders 
**From here we will refer to the buttons described in the [Introductory Guide](../GetStarted/NormalizingOneSpectrum.md).** 

### Useful Navigation Buttons:

The following buttons will be useful for completing this tutorial:

* To see the wavelength ranges of the orders, you can open the polynomial degree panel using the `set poly. degree...` button. When you hover the mouse over a specific order, a box will appear indicating the wavelength range of the order. Also, colors in the panel match the plot-colors of their respective orders. 
* Use the `Zin` (zoom in) and `Zout` (zoom out) buttons to adjust the plot so that the individual lines are visible. It is recommended that you also visualize the width of the entire spectral order you are working on.
* To view more of the depth of the line when they appear small, use the zoom button and select an enclosed area that is the width of the window but a much smaller height, with the height starting at the continuum and going down to the bottom of the spectral line. 
* If you have an ok horizontal range, but are too zoomed in or too zoomed out vertically, you can try the `auto-y` button, to automatically scale the vertical axis.
* Matplotlib's default panning and zooming buttons get overridden by the `include range` and `exclude range` functions. Make sure they are not on when you’re exploring the spectra. However, the custom buttons for zooming, panning, and auto-scaling will still work, as will the keyboard shortcuts.
* To see the outcome of any changes you’ve enacted, you have to click `fit cont`.

### Setting global fitting parameters:

Note that this tutorial uses the default 500 km s$^{-1}$ search bin width.

### Telluric/ Noise Dominated orders: 
In some instances, the observations may have have orders that have very low signal to noise ratios or which are dominated by telluric lines (sharps lines created from absorption by Earth's atmosphere). It can be difficult to try to apply the full set of fitting procedures on these orders to get them flat. Moreover, these orders contain very little useful information about the real spectrum of the star. 

In these cases, the priority should be to achieve a reasonably well-behaved continuum polynomial fit. You can try reducing the order of the polynomial degree by clicking the `set poly. degree...` button. For our example, change orders 1, 2, 38, 39, 40 to have smaller values (like 1 or 2). 

### Removing Spectral lines from the fit:

It is important to make sure that spectral lines are not contributing to the polynomial fit. None of the "good" continuum points should appear within spectral lines. Use the `exclude range` function to remove points within spectral lines from the fit. For example, let’s remove the $H \beta$ line at ~486 nm.

```{image} ../normplot_images/486_a.png
:alt: 486_a
:class: bg-primary mb-1
:width: 600px
:align: center
```
Select the region using the `exclude range` button: 

```{image} ../normplot_images/486_b.png
:alt: 486_a
:class: bg-primary mb-1
:width: 600px
:align: center
```
After pressing `fit cont`:

```{image} ../normplot_images/486_c.png
:alt: 486_a
:class: bg-primary mb-1
:width: 600px
:align: center
```
The wings of the $H \beta$ line should no longer be arching above the continuum since we exluded the problem points from the fitting routine!

If your exclude region is too big, you can use the `include range` button to bring useful portions of the spectrum back into the fit.

Now, try removing the spectral lines at the following wavelengths: 
383 nm, 397 nm, 410 nm, 434 nm ($H \gamma$), 438 nm, 646 nm, 656 nm ($H \alpha$), 850 nm, 860 nm, 866 nm, 875 nm, 810 nm, 901 nm, and 905 nm.

If we adjust the width of the search bins later, we will not have to worry about new points appearing on the same line. 

In some instances, we will have spectral lines that exist on the overlapping edge of two orders like $H \alpha$ and $H \beta$. For some lines, this can be managed by having the `fill order edge gaps` option turned on. Then, for that excluded region, the program will try to complete the fit using the nearest fit point from the neighboring order.


### Fit Points in Telluric Lines:
You may occasionally find some fit points being placed on telluric lines. This usually only happens in places where a lot of tellric lines blend together, so that there is no good continuum in a region of spectrum. We exclude those regions in an order, again with the `exclude range` button (which is especially useful if you choose to adjust the bin widths later on). For example let’s exclude the region around 759 nm. 

```{image} ../normplot_images/telluric_a.png
:alt: 486_a
:class: bg-primary mb-1
:width: 600px
:align: center
```
Select the region using the `exclude range` button: 

```{image} ../normplot_images/telluric_b.png
:alt: 486_a
:class: bg-primary mb-1
:width: 600px
:align: center
```
After pressing `fit cont`:

```{image} ../normplot_images/telluric_d.png
:alt: 486_a
:class: bg-primary mb-1
:width: 600px
:align: center
```
The wings of the telluric region around 759 nm are no longer arched and there are no fit points within the blended region.

### Fit Points on Noise: 
You may occasionally find a fit point resting in the noise level on the continuum. To correct this, adjust the **average length** and increase it to a slightly larger number. There is a small window of width that we specify that is moving across the spectrum one data point at a time and averaging the flux at each step about the center of that window. Adjusting the average length essentially smooths the spectrum contained in that bin, but does not directly affect the data. In effect, this will decrease the impact of the noise on where the fit point is placed. 

## Saving Normalized Spectra 

The chosen parameters (after fitting) can be saved to the files:  `exclude.dat`, `poly-deg.dat`, and `params.dat`.  These files can later be loaded to start from where you left off, or for a good initial guess if you are normalizing similar spectra.  Usually it's safest to always click this button before closing the main window, just in case you want to tweak a normalization later.  (Closing the main window always automatically saves the normalized spectrum.)

###  Output files
Once the main window is closed the program will save the normalized spectrum to `[observation_file].norm`.
```{note}
The output will use the normalization and parameters of the last update. It is advised that you use the `fit cont.` button before closing the window.
```

If you clicked the `save params`, the program will also write out the fitting information in the files, mentioned above. It's is useful to make a copy of these files with a name unique to the data. 

### Running normPlot with previous parameters
When running normPlot again, it will check for the prameter files and try to read them if they exist: `exclude.dat`, `poly-deg.dat`, and `params.dat`.  Alternatively, you can run normPlot using copies of those files with specified names from the command line:
```
normplot -e [exclude_file] -d [poly-deg_file] -c [params_file] [observation_file]
```
or within Python:
```
import normPlot
normPlot.normplot('[observation_file]',
                  excludeRegionName='[exclude_file]', 
                  polynomialsName='[poly-deg_file]', 
                  paramsName='[params_file]')
```

You can also take the parameter files and apply them to spectra without having the reopen the interface.  This uses the `-b` flag from the command line:
```
normplot -b -e [exclude_file] -d [poly-deg_file] -c [params_file] [observation_file]
```
or the `batchMode=True` argument in Python:
```
normPlot.normplot('[observation_file]',
                  excludeRegionName='[exclude_file]', 
                  polynomialsName='[poly-deg_file]', 
                  paramsName='[params_file]',
                  batchMode=True)
```
This useful for normalizing a large number of observations of the same object. 




