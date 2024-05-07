# Quickstart with the command line

This guide covers some of the more frequently used features of SpecpolFlow. It is aimed at people who are already familiar with spectropolarimetric data analysis and LSD. For more in-depth guides see [Normalizing one spectrum with the interactive tool](NormalizingOneSpectrum.md), [From normalized spectrum to Bz measurement](OneObservationFlow_Tutorial.ipynb).

To begin, install SpecpolFlow with pip, see [Installation](../About/Installation.md).

## Converting file formats

The SpecpolFlow tool-set uses the text .s file format for spectra. This usually consists of columns of wavelength, Stokes I, Stokes V, Null 1, Null 2, and sigma (wavelength, Stokes I, and sigma is also supported).  Converters for a few more common, instrument specific, .fits formats are provided.

For ESPaDOnS .fits spectra from the CADC use:
```
spf-fitstos-espadons observation.fits
```
This will generate a set of output files: the unnormalized spectrum (ending in u.s, e.g. `observationu.s`), automatic pipeline normalized spectrum (ending in n.s, e.g. `observationn.s`), and a header with information from the data reduction pipeline (ending in .out, e.g. `observation.out`).

## Normalizing spectra

Reliable continuum normalization is needed before running LSD. Several normalization tools can do this.  SpecpolFlow includes the optional tool normPlot, which can handle the .s format with polarimetric information. **Before using this, make sure you used the option to install normPlot!** (`pip install normPlot`)

Run the interactive tool on an input observation like:
```
normplot observationu.s
```
On closing the window, the normalized spectrum will be saved to a file `[observation_name].norm`. (It can be helpful to rename the normalized file like `mv observationu.s.norm  observationu-n.s`)

The tool normalized spectra by fitting low degree polynomials through continuum points, and operates on spectral orders individually. In the normalization Stokes I, V, and the nulls are all divided by the continuum polynomial. There are tooltips in the program to briefly explain the controls.

The tool can also run in a non-interactive 'batch' mode, just reading control parameters from files like:
```
normplot observationu.s -b -e exclude.dat -d poly-deg.dat -c params.dat
```

A more detailed guide is [Normalizing one spectrum with the interactive tool](NormalizingOneSpectrum.md).  For details on the command line parameters run `normplot -h`

## Making and cleaning a mask

### Making the mask

To make a line mask, start from a VALD 'extract stellar' line list in 'long' format.  The line data, and depths, are taken from the VALD list.  To convert the list to a mask run:
```
spf-makemask line_list.dat line_mask.dat
```

The line list can also be trimmed to only lines deeper than some depth, and in some wavelength range, like:
```
spf-makemask -d 0.2 -w1 450 -w2 650 line_list.dat line_mask.dat
```
For a full list of options run `spf-makemask -h`

### Cleaning the mask

It is strongly recommended to clean the line mask before using it: removing lines that are blended with telluric features, or big broad spectral features like Balmer lines that LSD can't reproduce.  SpecpolFlow has an interactive graphical tool for this, which can be run like:
```
spf-cleanmask line_mask.dat observationu-n.s line_mask_clean.dat
```
This displays the line mask, a reference observation, and the LSD model's fit to the observation. Lines can be selected to remove from the mask (red) or include (blue).  The tool also allows for fitting line depths, assuming the current LSD profile is correct.  However, line depth fitting should be treated with some caution since there is an intrinsic degeneracy between LSD profile amplitude and line mask depths.

The mask cleaning can also be run in a non-interactive batch mode, using a file of exclude regions, like:
```
spf-cleanmask -b -e exclude_regions.dat line_mask.dat observationu-n.s line_mask_clean.dat
```
For a full list of command line parameters run `spf-cleanmask -h`.  For more details on using the program see [How to clean masks with the interactive tool](../Tutorials/3b-MaskUI_Tutorial.md).


## Running LSD

The LSD code in SpecpolFlow, LSDpy, can be run like:
```
lsdpy -m line_mask_clean.dat observationu-n.s lsd_profile.dat
```
LSDpy will try to read a set of input parameters from a file `inlsd.dat`, and output some diagnostic information to the terminal.  For a full list of command line parameters run `lsdpy -h`.

If the `inlsd.dat` file does not yet exist, LSDpy will generate a template file and halt.  You can then edit the template for parameters specifically for your observations, and re-run the code.

In the `inlsd.dat` file, pay particular attention to the  "start and end velocity", "pixel size in velocity", and "mask/profile normalization parameters" values.  It can also be helpful to set "Save LSD model spectrum?" to "`1  outModelSpec.dat`", to generate the LSD best fit model spectrum for direct comparison with an observation.

The information printed to the terminal can be useful, and includes information about the line mask, information about the treatment of error bars, and information about detection statistics in Stokes V and the null.

## Processing LSD profiles

### Plotting LSD profiles
For quickly plotting LSD profiles:
```
spf-plotlsd profile1.dat
```
Or plotting multiple profiles with a legend displaying file names:
```
spf-plotlsd -l profile1.dat profile2.dat profile2.dat
```

### Calculating Bz
For calculating longitudinal magnetic fields (Bz), and detection statistics, for Stokes V and the null, you can use spf-bz like:
```
spf-bz -v -20 +30 -g 1.2 -l 500 profile1.dat
```
This example uses a velocity integration range from -20 km/s to +30 km/s, an effective Landé factor of 1.2, and a wavelength of 500 nm, to calculate Bz from profile1.dat.  The False Alarm Probability (FAP) inside the -20 to +30 km/s line is also calculated, and printed to terminal.  

To correctly choose the velocity range to use, it can be helpful to plot the LSD profile and integration range with the `-p` option:
```
spf-bz -p -v -20 +30 -g 1.2 -l 500 profile1.dat
```

In general the effective Landé factor and wavelength used here should be the normalizing values used for the LSD calculation.  Within the SpecpolFlow tools, LSDpy saves those values to the header of the LSD profile, and `spf-bz` can read them, so you have the option of just using:
```
spf-bz -v -20 +30 profile1.dat
```
(This won't work for LSD profiles from other codes).
Multiple files can be processed at once:
```
spf-bz -v -20 +30 profile1.dat profile2.dat profile3.dat
```
For a full list of command line parameters run `spf-bz -h`.  For a detailed explaination of these calculations see [Calculating the longitudinal field and False Alarm Probability](../Tutorials/6-CalculateBz_Tutorial.ipynb).


### Calculating radial velocity

There is a simple tool for quickly calculating radial velocities from LSD profiles:
```
spf-rvfit -v -100 +100 profile1.dat
```
This fits a Gaussian to the LSD profile (in this example using the -100 to +100 km/s range), and takes the centroid as the radial velocity.  The fit to the LSD profile can be plotted with the `-p` option.  Multiple files can be processed at once:
```
spf-rvfit -v -100 +100 profile1.dat profile2.dat profile3.dat
```
Command line parameters can be printed with `spf-rvfit -h`.  While this Gaussian fitting approach is convenient and often reliable, alternative approaches may be more optimal for your data.
