# Normalizing a Single Spectrum with NormPlot

## About NormPlot2.4 

This current version of the interactive UI is written in Python 3.8. Its major dependecies are **tkinter** and **matplotlib** (if they are not installed python will issue a running error). NormPlot can be utilized as a preliminary step in the `SpecpolFlow` pipeline to normalize the continuum of 1D stellar spectra (post-reduction). It is designed to work well with echelle spectra. This interactive UI was built to allow you to work on one sepctral order at a time in which you essentially establish "good" continuum points and fit a low order polynomial through them.


## Using the GUI

```{note}
 Prior to using the normalization GUI, the file containing the Stokes data should be converted from fits to .s format that follows the LibreESPRIT format. Please see the tutotial or documentation regarding the use codes for conveting to .s files as some have been provided for commonly used spectropolarimeters.  
```


Assuming you have installed the program via pip, the interactive normalization UI can be started from the terminal with the command 

```
normplot [observation_file]
```

`````{admonition} Command line alternative (if NormPlot is not pip installed)
:class: tip 
Use the command  
python normPlot2.py [observation_file]
`````

Once executed, you should see the following window appear

```{image} ../normplot_images/user_guide_gui.png
:alt: gui_snapshot
:class: bg-primary mb-1
:width: 800px
:align: center
```

The top panel of this window displays the spectral orders and polynomial fits for the black scatters points on the spectrum of each order. The black fitting points are computed by dividing an order into consecutive search bins of equal width and extracting the mean flux within each bin. The panel below is a preview of the normalization product that would be output by the code. The axes are synced and allow the user to comparatively inspect  a spectral order and the respective normalized spectrum within the same range of wavelength. The lower panel updates anytime the new fit is explicitly computed. **The wavelength  units are in nanomneters.**

### Button Functions in the Interactive UI

There are controls for navigating the GUI plot:
`auto Scales`, `auto-y`, `Zin` (zoom in), `Zout` (zoom out), and panning buttons (left, right, up, and down).The python default below these buttons also work for navigating the frame home position, revert step, redo step, pan, zoom, and snapshot.

There are different tools for deriving an optimal fit of the individual orders. Here, we explain the other functions and features that can be utilized in the fitting process. 

**Fit cont.**: 
Clicking this button will apply any parameter changes made and update the continuum fit. This will simultaneously update the final normalized spectrum. 


**set output params**:  

```{image} ../normplot_images/user_guide_set_output.png
:alt: gui_set_output
:class: bg-primary mb-1
:width: 400px
:align: center
```

Here you can determine the structure of the spectrum in the output file. If you select **merge spectral orders**, the final spectrum will be continuous rather than separate orders (if separate there will be wavelength overlap between consecutive orders). In the process of merging the orders, the program will choose the spectrum of the upper order. The lower panel in Figure 1 will change accordingly, so you will be able to preview the merged normalized spectrum. The wavelength units can also be scaled. If converted to angstroms (set the scale to 10), the wavelengths can be converted between air and vacuum wavelengths. 

**set poly. degree**:

```{image} ../normplot_images/user_guide_polyfit_params.png
:alt: gui_set_output
:class: bg-primary mb-1
:width: 200px
:align: center
```

Here the degree of the polynomial fits (of the black bin points)  can be set for each individual order. The color of the text corresponds to the color by which that order is plotted. If you hover your mouse over one of the inserts, text will pop up indicating the wavelength range for that given order. This will help in searching for a specific order on the plot. After clicking apply and then `fit cont.`, you will see a change in the polynomial plotted over the order and in the normalized spectrum in the region corresponding to that order. 

**fill edge order gaps**: 
This check box enables a feature such that if a spectral order ends in an exclude region, the algorithm can take a fit point in the adjacent spectral order and use that to constrain the continuum polynomial

**average length**:

Within all the search bins, you can apply boxcar averaging to smoothen the data before extracting the point (mean flux). The input sets the width (ie. number of neighboring data points on the spectrum) that will be used in the averaging. This is most advantageous to correct for fit points that favor noise. This smoothening is only in effect when computing points for the fits of the orders, so you will not produce a smoothed normalized spectrum.     

**srch. bin (km/s)**:

This value sets the width of the search bins. Recall, that there is only one point extracted per bin. As you increase the width of the bins, you decrease the number of fit points for each order and vice versa. 

**exclude range**:

You can select specific ranges of wavelengths to be excluded from the fitting process. After clicking this option and going into either panel, the first click will set the starting edge of the exclusion range. A box made of dashed lines will appear with one edge fixed on where the first click. The second vertical edge will follow the mouse in the panel, and  all wavelengths within this box will be excluded. A second click will finalize this selection and after pressing fit. cont. you will see that enclosed portion of the spectrum change its color to black (this signifies that it is not included in the fitting). Any fit points in that enclosed region will also disappear. 

**include range** 
This provided a way to reverse the exclusion function. You can selected ranges with in excluded regions to be unexcluded when the normalization is updated. 

**hide norm**
You have the option to remove the lower panel (or return it) at your leisure. The upper panel will expand to the width of the full window. This feature can be useful for those doing the normalization of small computer screens, and would rather have fuller viewing of the original spectrum and fits. This is convenient for users operating the program on small screens.

**save params**: 
The chosen parameters (post fitting) can be saved to the files:  **exclude.dat**, **poly-deg.dat**, and **params.dat**.

###  Outputs 
Once the window is closed the program will save the normalzied spectrum to **observation_file.norm**.
```{note}
The output will keep the normalzition and paramters of the last update. It is advised that you use the **fit cont.** function before closing the window.
```

The program will also write out the fitting information in the files mentioned in the **save params** button description. It's is advisable to make a copy of these files with a name unique to the data. 

You can then take the updated parameter file and apply them to spectra without having the reopen the interface by executing
```
python normPlot2.py observation_file.s -b
```
This useful for normailizing multiple observations of the same object. 



