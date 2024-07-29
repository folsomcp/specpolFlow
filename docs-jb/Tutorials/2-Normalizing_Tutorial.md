# How to Normalize echelle spectra with NormPlot

:::{Warning}
Still under construction. Some commands and other details may be out of date.
:::

## Introduction

This tutorial introduces the procedure of normalizing echelle spectra and the use of NormPlot. The interactive normalization GUI is a part of the SpecpolFlow pipeline as a step that proceeds the polarimetric analysis. In this introduction to spectra normalization, we will be using ESPADonS observations of the star $\xi^1$ CMa. This is the same observation that is used in the LSD and Mean Longitudinal Field tutorials, so skipping this tutorial will not prevent you from advancing further in the pipeline (If you'd like, you're welcome to use normalized spectrum you create while doing this tutorial and see how the magnetic field measurements compare). 

::::{margin}
:::{note}
This current version of the tutorial focuses on executing NormPlot from the command line. Refer to the user guide for executing the same commands from a Jupyter Notebook or python script. 
:::
::::


## Converting from Fits to text
Let's proceed by converting the format for the data files from '.fits' to the '.s' text format. The data file needs to be converted in order to be used in NormPlot. These processed data files (which can be found in CADC archives) are formatted as follows: `i.fits` contains spectroscopic data (usually 4 spectra for each polarimetric observation), and `p.fits` contains the spectopolarimetric data. We want to convert the `p.fits` files and use the Stokes I spectrum in the normalization code. For this tutorial, we will use the ESPaDOns observation, `2378200p.fits`.

::::{margin}
:::{note}
ESPaDOns `p.fits` files contain a main data table with 24 columns. Here we are mostly interested in columns 19-24, which contain the unnormalized spectrum (without a telluric velocity correction), and have columns for wavelength, intensity, polarization, null 1, null 2, and errorbars.
:::
::::

In the command line, write the command: 

 ``` 
 spf-fitstos-espadons inputfile-p.fits
 ```

but replace the `inputfile-p.fits` with the name of our ESPaDOnS file `2378200p.fits`. The script will return `2378200pn.s` (normalized data) and `2378200pu.s` (unnormalized data) in the same directory as the '.fits' file. 

This example is specific to ESPaDOnS observations from the standard Upena/LibreESPRIT pipeline at CFHT. For observations from other instruments or other data reduction pipelines, you will need a different script to convert FITS files to the `.s` format. In some cases you may need to write your own conversion script; converters for a few instruments are discussed in [this Tutorial](1-ConvertToSFiles_Tutorial.ipynb), and documented in the [API](https://folsomcp.github.io/specpolFlow/API/Converters_API.html).


## Normalization GUI 

**Our goal in this process is to fit a curve through the star’s continuum in every spectral order. With that, we can normalize the continuum by dividing the observation by the fit curve. This creates a ‘common line’ (a flat horizontal line at y = 1) about which we can consistently measure the properties of spectral lines.** 

Now let’s use the normalization GUI on our data. This normalization code utilizes the text format spectrum that would have been produced in the first steps (it was built around LibreESPRIT reduced data). You can run the code in the command line with: 

```
normplot <observation_file.s>
```

Once executed, a window will pop up and you will see the spectrum separated into orders, indicated by different colors. There are polynomial curves plotted, which are fits to the continuum for each order, colored differently for each order. There are also black points plotted, which are the points the used for the polynomial fit.

```{image} ../normplot_images/gui_open.png
:alt: gui_open
:class: bg-primary mb-1
:width: 600px
:align: center
```

In the bottom left corner, there is a toolbar with commands from matplotlib (panning, zooming, returning to default, saving a figure). Above that are a set of buttons with some custom commands for automatically scaling both axes, auto-scaling just the y-axis, zooming, and panning. There are keyboard shortcuts for these commands (arrow keys for panning, 'a' for auto-scaling, 'A' for auto-scaling the y-axis, 'z' for zooming on a selected region, 'i' and 'o' for zooming in and out).

**The fit points**: The idea here is to select a relatively small number of 'good' continuum points, and use those to fit a polynomial describing the continuum of the observation. These 'good' points are shown as black dots on the plot. To find those 'good' points, each order is broken down into a set of consecutive search bins with respect to wavelength. The bins are equal in width (in velocity). In each bin, the highest flux point is selected as being most likely in the continuum (i.e. not in an absorption line). However, to deal with noise in the observation, a running average is applied to the spectrum before selecting the highest point in the search bin. This averaging helps ensure that the selected highest point in a bin is not highest because of a peak in the noise, but instead highest because of real features in the spectrum.
Collectively, the black dots within an order are the points that are used to fit a polynomial function for the continuum of that order. The goal of this process is selecting black dots representative of the continuum flux so that when we normalize, the continuum is near-uniformly 1. In other words, you want a flat line at f = 1 with various spectra lines added onto it. 


## Normalizing the Individual Orders 
**From here we will refer to the buttons described in the [Introductory Guide](../GetStarted/NormalizingOneSpectrum.md). This guide addresses items that need attention in the normalization process.** 

### Navigation items to remember:

A few things that will be useful when doing this tutorial:

* If you want to see the wavelength ranges of orders, you can open the polynomial degree panel (with the `set poly. degree...` button). When you hover the mouse over a specific order, a box will appear indicating the wavelength range of the order. Also, colors in the panel match the plot-colors of their respective orders. 
* It is useful to zoom in enough that you can see individual lines, but keep the plot wide enough that you can see the width of the entire order you are working on.
* To view more of the depth of the line when they appear small, use the zoom button and select an enclosed area that is the width of the window but a much smaller height, with the height starting at the continuum and going down to the bottom of the spectral line. 
* If you have an ok horizontal range, but are too zoomed in or too zoomed out vertically, you can try the `auto-y` button, to automatically scale the vertical axis.
* Matplotlib's default panning and zooming buttons get overridden by the `include range` and `exclude range` functions. Make sure they are not on when you’re exploring the spectra. However, the custom buttons for zooming, panning, and auto-scaling will still work, as will the keyboard shortcuts.
* To see the outcome of any changes you’ve enacted, you have to click `fit cont.`.

### Setting global fitting parameters:

Say something about the search bin width and average length for starting values?

**For now we will keep the default 500 km s$^{-1}$ search bin width.** 

### Telluric/ Noise Dominated orders: 
Let’s get through some of the less rigorous things first. In some instances, your observations may have have orders that have very low signal to noise ratios or are dominated by telluric lines (sharps lines created from absorption by the atmosphere). It can be difficult to try to apply the full set of fitting procedures on these orders to get them flat. Moreover, these orders contain very little useful information about the real spectrum of the star. 

In these cases, where you aren't planning to use the order anyway, we mostly just want the continuum polynomial to be reasonably well behaved. One of the most direct fixes to this is to reduce the order of the polynomial degree. In our example, click `set poly. degree...` button, then change orders 1, 2, 38, 39, 40 to have lower values. 1 or 2 might be a good value. 


### Removing Spectral lines from the fit:
Remember, we essentially are making a common line that the spectral lines extend down from, so we must do our best to make sure they aren’t contributing to the fitting. This means we don’t want any of our fit points resting on spectral lines. Wherever we see fit point on spectral lines, we use the `exclude range` function to remove those line from being options to place fit points. For example, let’s remove this line at ~383 nm.
If your exclude region is too big, you can use `include range` to regions to bring back portions of the spectrum that you would like back.
Now, there are a few other lines that we need to tackle at: 397, 410, 434 ($H \gamma$) , 438, 486 ($H \beta$), 646, 656($H \alpha$), 850, 860, 866, 875, 810, 901, and 905nm.

**figure?**

One advantage of this is that if we adjust the width of the search bins later, then you will not have to worry about new points appearing on the same line. 
In some instances, you will have spectral lines that exist on the overlapping edge of two orders like H_alpha and H_beta. For some lines, this can be managed by having the `fill order edge gaps` on. Then for that excluded region, the program will try to complete the fit using the nearest fit point from the neighboring order.


### Fit Points in Telluric Lines:
You may occasionally find some fit points being placed on telluric lines. This usually only happens in places where a lot of tellric lines blend together, so that there is no good continuum in a region of spectrum. We exclude those regions in an order, again with `exclude range` (especially useful if you choose to adjust the bin widths later on). For example let’s exclude the region about 759 nm. 

**figure?**

### Fit Points on Noise: 
You may occasionally find a fit point resting in the noise level on the continuum. To try to correct this, we can adjust the **average length** and increase it to a slightly larger number. What’s happening here is that there’s a small window of width that we specify that is moving across the spectrum one data point at a time and averaging the flux at each step about the center of that window. This essentially smoothens the spectrum contained in that bin, but this does not directly affect the actual data we use. It will lessens the impact of the noise on where the fit point is placed. 


## Saving Normalized Spectra 

You can use merge the product so it keeps the spectrum of the upper order, or you can keep it unmerged and use your own discretion in how you use the overlapping regions. 

Finally save in output and make a copy of the file should you need to revisit the normalization later. 

You can also use saved input parameters from the command line. Run the code with `python normPlot2.py -h` for some extra information about that. 



