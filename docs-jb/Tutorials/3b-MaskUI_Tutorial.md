# How to clean masks with the interactive tool

## Introduction

This tutorial will introduce you to the interactive mask cleaning tool called `cleanMaskUI`. This tool allows regions of lines to be easily excluded and included in the LSD line mask. This is particularly useful in the case of stars with a lot of emission: since emission lines have a different shape than the typical absorption lines, they should not be included in the line averaging since they will distort the resulting LSD profile. For this tutorial we will use the star HD 164284, a Classical Be star with a spectrum that is contaminated by emission lines from a disk. 

`cleanMaskUI` can be run from either a Python file, a cell within a Jupyter Notebook, or from the command line.

From a cell we can call the `cleanMaskUI` function. This function takes in:
- The file name of your line mask, or the corresponding Mask object
- The file name of an observed spectrum `.s` file, for comparison and LSD calculations
- The file name for saving the new cleaned mask (optional; if not specified it will use the input mask name with '.clean' added)
- The exclude mask region file, which contains the set of regions to be excluded from the mask (optional)

For running the tool from the terminal use:
```
spf-cleanmask input_mask_file reference_observation.s cleaned_mask_file
```

For running the tool from inside a Python or Jupyter Notebook use:
```
import specpolFlow as pol
pol.cleanMaskUI('input_mask_file', 'reference_observation.s', 'cleaned_mask_file')
```

Running the code opens a new window, with a GUI that will look like this:
```{image} ../MaskUI_images/openUI.png
:alt: openUI
:class: bg-primary mb-1
:width: 600px
:align: center
```

The vertical lines indicate the locations and depths of spectral lines in the input line mask. Lines colored blue indicate that they are being used (with `iuse=1`), and lines colored red indicate that they are not being used (with `iuse=0`). The black line is the observed spectrum, and the purple line is the model LSD spectrum (i.e. the spectrum created by convolving the LSD profile with the line mask). 

## Navigating

The GUI's navigation controls are located in the bottom left. Hovering the mouse over the different buttons will provide a brief description of what that button does. To zoom in on a specific region we can select the `zoom` button, or the magnifying glass symbol. Now we can click and drag on the plot to create a dashed box around a region. When we release the mouse, the plot will rescale to the selected region.

To return to the original view we can click `auto.` button, or the home button. To only scale the vertical axis use the `auto-y` button. We can also click the back arrow to go back to the previous view. 

There are keyboard shortcuts for many of these buttons (arrow keys pan, `a` auto-scales, `i` and `o` zoom in and out).  Mouse over a button to see the keyboard shortcut.

The figure below shows a zoom in of a region including a line with emission present. 

```{image} ../MaskUI_images/zoomIn.png
:alt: zoomIn
:class: bg-primary mb-1
:width: 600px
:align: center
```
We can also click the `Zin` and `Zout` buttons to quickly zoom in or out by a small amount. 

## Including and Excluding Lines

Lets now manually exclude a region around that emission line. On the bottom right click the `exclude lines` button. The button should now have a dotted box around it indicating that it is active, as shown below. 

```{image} ../MaskUI_images/excludelinesButton.png
:alt: excludelinesButton
:class: bg-primary mb-1
:width: 300px
:align: center
```

We can now click on the plot to create a vertical dotted line indicating one edge of the region. We can then click on another part of the plot to finish selecting the region. The region selection is shown below.

```{image} ../MaskUI_images/excludingRegions1.png
:alt: excludingRegions1
:class: bg-primary mb-1
:width: 600px
:align: center
```

Once we click a second time to close the region, you should see that all spectral lines within the selected region now turn red. (You can also right click during the selection, after the first click, to cancel.)

```{image} ../MaskUI_images/excludingRegions2.png
:alt: excludingRegions2
:class: bg-primary mb-1
:width: 600px
:align: center
```

Once you have finished selecting regions to exclude, don't forget to click the `exclude lines` buttons again to deselect that mode. You should see the dotted box around the button disappear.  Note that the custom buttons for panning and zooming will work while the 'exclude lines' mode is active, but the default matplotlib toolbar at the very bottom does not.

The above procedure works exactly the same for the `include lines` button. That button will add lines back into the mask, in other words it will turn the red (excluded) lines back to blue (included).  By default, lines in some commonly problematic regions (telluric or Balamer line regions) are excluded, but in some cases it may be helpfull to include those lines.

## Updating and Plotting the LSD Profile

Once we have selected regions to exclude we can update the test LSD profile by clicking the `update LSD` button in the bottom right. This will update the purple model spectrum line as shown below and update the output cleaned mask. 

```{image} ../MaskUI_images/updateLSD1.png
:alt: updateLSD1
:class: bg-primary mb-1
:width: 600px
:align: center
```

We can also adjust the input parameters for the test LSD profile by selecting the `LSD param.` button in the bottom right. This will open another window titled `Set LSD parameters` as shown. 

```{image} ../MaskUI_images/updateLSD2.png
:alt: updateLSD2
:class: bg-primary mb-1
:width: 200px
:align: center
```

For this example, I will change the starting and stopping velocity to $\pm$ 800 km/s and I will select the `remove closely spaced lines` and `plot profile` boxes. The latter box will actually plot the test LSD every time we press the `update LSD` button. (In a Jupyter Notebook this LSD plot should be visible back in the cell output where we first ran the `cleanMaskUI` function.)

```{image} ../MaskUI_images/updateLSD3.png
:alt: updateLSD3
:class: bg-primary mb-1
:width: 200px
:align: center
```

We can now close the `Set LSD parameters` popup window and select `update LSD` to update the model spectrum and create a LSD plot. 

```{image} ../MaskUI_images/updateLSD4.png
:alt: updateLSD4
:class: bg-primary mb-1
:width: 400px
:align: center
```

Since we only removed one emission line so far, and did not change any other lines, the LSD profile looks quite bad. However, after doing a more thorough pass we can obtain an LSD profile that is much more usable. 

```{image} ../MaskUI_images/updateLSD5.png
:alt: updateLSD5
:class: bg-primary mb-1
:width: 400px
:align: center
```

:::{Warning}
Still under construction.  More details forthcoming...
:::

<!-- Fit depths (TODO) -->

