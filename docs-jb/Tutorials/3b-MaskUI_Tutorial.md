# How to clean masks with the interactive tool

## Introduction

This tutorial will introduce you to the interactive mask cleaning tool called `cleanMaskUI`. This tool allows regions of lines to be easily excluded and included in the LSD line mask. This is particularly useful in the case of stars with a lot of emission: since emission lines have a different shape than the typical absorption lines, they should not be included in the line averaging since they will distort the resulting LSD profile. For this tutorial we will use the star HD 164284, a Classical Be star with a spectrum that is contaminated by emission lines from a disk. 

`cleanMaskUI` can be run from either from the command line, or from a Python file, or a Jupyter Notebook.

The `cleanMaskUI` function takes:
- The file name of your line mask, or the corresponding Mask object
- The file name of an observed spectrum '.s' file, for comparison and LSD calculations
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
```{image} MaskUI_images/openUI.png
:alt: openUI
:class: bg-primary mb-1
:width: 600px
:align: center
```

The vertical lines indicate the locations and depths of spectral lines in the input line mask. Lines colored blue indicate that they are being used (with `iuse=1`), and lines colored red indicate that they are not being used (with `iuse=0`). The black line is the observed spectrum, and the purple line is the model LSD spectrum (i.e. the spectrum created by convolving the LSD profile with the line mask). 

## Navigating

The GUI's navigation controls are located in the bottom left. Hovering the mouse over the different buttons will provide a brief description of that button's action. To zoom in on a specific region, select the `zoom` button. 
A region can be selected by clicking on the plot, moving the mouse to create a box around the desired region. Then clicking the mouse a second time will rescale the plot to the selected region.
Or you can use the the click magnifying glass symbol for the matplotlib zoom tool (click and drag to select a region, release the mouse button will rescale to that region).

To return to the original view, click the `auto` button (or the matplotlib home button). To only scale the vertical axis, use the `auto-y` button. Clicking the back arrow will return to the previous view. 

There are keyboard shortcuts for many of these buttons (arrow keys to pan, `a` to auto-scale, `i` and `o` to zoom in and out). Hover the mouse over a button to see the associated keyboard shortcut.

The figure below shows a zoomed-in view of a region including a line with emission present. 

```{image} MaskUI_images/zoomIn.png
:alt: zoomIn
:class: bg-primary mb-1
:width: 600px
:align: center
```
We can also click the `Zin` and `Zout` buttons to quickly zoom in or out by a small amount. 

## Including and Excluding Lines

Lets now manually exclude a region around that emission line. On the bottom right, click the `exclude lines` button. The button should now have a dotted box around it indicating that it is active, as shown below. 

```{image} MaskUI_images/excludelinesButton.png
:alt: excludelinesButton
:class: bg-primary mb-1
:width: 300px
:align: center
```

Click on the plot to create a vertical dotted line indicating one edge of the region. Then, click on another part of the plot to finish selecting the region. The region selection is shown below.

```{image} MaskUI_images/excludingRegions1.png
:alt: excludingRegions1
:class: bg-primary mb-1
:width: 600px
:align: center
```

After the second click to close the region, all spectral lines within the selected region should now turn red. (You can also right click during the selection, after the first click, to cancel.)

```{image} MaskUI_images/excludingRegions2.png
:alt: excludingRegions2
:class: bg-primary mb-1
:width: 600px
:align: center
```

Once you have finished selecting regions to exclude, don't forget to click the `exclude lines` button again to deselect that mode. You should see the dotted box around the button (and lighter color) disappear.  Note that the custom buttons for panning and zooming will work while the 'exclude lines' mode is active, but the default matplotlib toolbar at the very bottom does not.

The above procedure works exactly the same for the `include lines` button. That button will add lines back into the mask, in other words it will turn the red (excluded) lines back to blue (included).  

By default, lines in some commonly problematic regions (around telluric lines and Balamer lines) are excluded.  In some cases it may be helpful to include those lines. In other cases you may want exclude additional wider regions.

## Updating and Plotting the LSD Profile

Once we have selected the regions to exclude, we can update the test LSD profile by clicking the `update LSD` button in the bottom right. This will update the purple model spectrum line as shown below and update the output cleaned mask. 

```{image} MaskUI_images/updateLSD1.png
:alt: updateLSD1
:class: bg-primary mb-1
:width: 600px
:align: center
```

We can also adjust the input parameters for the test LSD profile by selecting the `LSD param` button in the bottom right. This will open another window titled `Set LSD parameters` as shown. 

```{image} MaskUI_images/updateLSD2.png
:alt: updateLSD2
:class: bg-primary mb-1
:width: 200px
:align: center
```

::::{margin}
:::{note}
In the `Set LSD parameters` window, the `remove closely spaced lines` option only removes lines inside the LSD calculation. It does not remove them from the mask being cleaned!
:::
::::
In this example we will change the starting and stopping velocity to $\pm$ 800 km s$^{-1}$, since this Be star has an extremely high $v\sin i$.  We also select the `remove closely spaced lines` and `plot profile` boxes. The latter box will plot the test LSD profile every time the `update LSD` button is pressed. In a Jupyter Notebook this LSD plot should be visible back in the cell output where the `cleanMaskUI` function was executed.

```{image} MaskUI_images/updateLSD3.png
:alt: updateLSD3
:class: bg-primary mb-1
:width: 200px
:align: center
```

Close the `Set LSD parameters` popup window and select `update LSD` to update the model spectrum and create an LSD plot. 

```{image} MaskUI_images/updateLSD-clean1.png
:alt: updateLSD-clean1
:class: bg-primary mb-1
:width: 400px
:align: center
```

Since only one emission line was removed, and since no other lines were changed, the LSD profile is still distorted, and noisy with large errorbars. However, after doing a thorough pass through the whole spectrum, a useful LSD profile can be obtained. 

Since the lines are relatively shallow in this star, we remove more lines blended with weaker telluric lines.

```{image} MaskUI_images/updateLSD-clean-telluric.png
:alt: updateLSD-clean-telluric
:class: bg-primary mb-1
:width: 600px
:align: center
```

We also remove some lines where the line depths are wrong and the lines aren't clearly seen in the observation.  In this star, that is partly due to the very high $v\sin i$.

```{image} MaskUI_images/updateLSD-clean-weak-nd.png
:alt: updateLSD-clean-weak-nd
:class: bg-primary mb-1
:width: 600px
:align: center
```

After going through the entire spectrum, removing more problem regions, we get a much better LSD profile.

```{image} MaskUI_images/updateLSD-clean-full.png
:alt: LSD profile for the fully cleaned mask
:class: bg-primary mb-1
:width: 400px
:align: center
```

## Automatically Adjusting Line Depths 

Sometimes it can be helpful to use empirical estimates of line depths for some lines in the mask.  The theoretical line depths could be wrong due to NLTE effects, errors in the atomic line data, or denaturation/broadening effects that weren't accounted for.  

In cases where many lines in the mask have depth problems, it is usually better to improve the theoretical depths going into the mask.  It is important to have the correct effective temperature and $\log g$ for the mask.  If the star has chemical peculiarities it is important to account for those too.  This can be done within a VALD extract stellar request.  You may also wish to consider more elaborate spectrum synthesis calculations for the line depths, particularly if NLTE effects are very important.

In generally empirical line depth estimates should be used with caution. There is a good reason for your choice of theoretical line depths (or there should be!).  Furthermore, there is an intrinsic degeneracy between the depths of the lines in the mask and the amplitude of the LSD profile.  In `cleanMaskUI` this is resolved by simply assuming the current LSD profile is an adequate approximation when trying to fit line depths.  This means that if your initial LSD profile (from the initial line mask) is too poor quality the depth fitting will produce poor results.  This also means that one should only fit carefully selected problem lines, not all the lines in the mask.

In our example Be star, fitting depths of some lines can be quite helpful.  Since $v\sin i$ is quite large there are relatively few well detected lines in the spectrum, so depth errors in some of those lines can have a large impact.

