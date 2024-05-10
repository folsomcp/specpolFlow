# How to clean masks with the interactive tool

## Introduction

This tutorial will introduce you to the interactive mask cleaning tool called `cleanMaskUI`. This tool is optional and allows regions of lines to be easily excluded and included. This is particularly useful in the case of stars with lots of emission. Emission lines have a different shape than the typical absorption lines and thus should not be included in the line averaging as the resulting LSD profile will be distorted. For this tutorial we will use the star HD 164284 to demonstrate the capabilities of this tool. HD 164284 is a Classical Be star and as such its spectrum is contaminated by emission lines from a disk. 

`cleanMaskUI` can be run from either a cell within a `.ipynb` file or from the command line.

From a cell we can call the `cleanMaskUI` function. This function takes in:
- The file name of your line mask, or the corresponding Mask object.
- The file name of the spectrum `.s` file
- The file name to save the new cleaned mask to (optional)
- The exclude mask region file which contains the set of regions to be excluded from the mask (optional)

Running the following code will open the GUI in a new window.
```
import specpolFlow as pol
pol.cleanMaskUI(mask_file,spectrum_file,cleaned_file)
```

Upon opening the GUI will look like this. 
```{image} ../MaskUI_images/openUI.png
:alt: openUI
:class: bg-primary mb-1
:width: 600px
:align: center
```

The vertical lines indicate the locations and depths of spectral lines in the input line list. Lines colored blue indicate that `iuse=1` so the line is being used, and red lines indicate `iuse=0` so the line is not being used. The blue line is the spectrum, and the purple line is the model spectrum (the spectrum created by convolving the LSD profile with the line mask). 

## Navigating

The GUI's navigation controls are located in the bottom left. Hovering the mouse over the different buttons will provide a brief description of what that button does. To zoom in on a specific region we can select the magnifying glass symbol. Now we can click and drag on the plot to create a dashed box around a region. When we release the mouse, the plot will rescale to the selected region.

To return to the original view we can click the home button. We can also click the back arrow to go back to the previous view. 

The figure below shows a zoom in of a region including a line with emission present. 

```{image} ../MaskUI_images/zoomIn.png
:alt: zoomIn
:class: bg-primary mb-1
:width: 600px
:align: center
```

We can also click the `Zin` and `Zout` buttons to quickly zoom in or out by a small amount. 

## Including and Excluding Lines

Lets now manually exclude a region around that emission line. On the bottom right click the `exclude lines` button. The button should now have a dotted box around it as shown below. 

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

Once we select a second line to close the region, you should see that all spectral lines within the selected region now turn red. 

```{image} ../MaskUI_images/excludingRegions2.png
:alt: excludingRegions2
:class: bg-primary mb-1
:width: 600px
:align: center
```

Once you have finished selecting regions to exclude in the current view, you MUST click the `exclude lines` buttons again to deselect that mode. You should see the dotted box around the button disappear. Now we can navigate to the next region we wish to exclude. 

The above procedure works exactly the same for the `include lines` button. That button will add lines back into the mask, in other words it will turn the red (excluded) lines back to blue (included).

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

For this example, I will change the starting and stopping velocity to $\pm$ 800 km/s and I will select the `remove closely spaced lines` and `plot profile` boxes. The latter box will actually plot the test LSD every time we press the `update LSD` button. This LSD plot should be visible back in the cell output where we first ran the `cleanMaskUI` function. 

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

Since we only removed one emission line and did not change any other lines, the LSD profile looks quite bad. However, after doing a more thorough pass we can obtain an LSD profile that is much more usable. 

```{image} ../MaskUI_images/updateLSD5.png
:alt: updateLSD5
:class: bg-primary mb-1
:width: 400px
:align: center
```

<!-- Fit depths (TODO) -->

