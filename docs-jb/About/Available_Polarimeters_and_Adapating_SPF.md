# General Comments on Compatibilities

SpecpolFlow is generally compatible with data in the [LibreESPRIT](https://www.cfht.hawaii.edu/Instruments/Spectroscopy/Espadons/Espadons_esprit.html) format. 
This is usually a text file with columns for wavelength, intensity, polarization, two polarimetric nulls, and uncertainties.

For spectra in a text file with columns for wavelength, intensity, and uncertainties, the polarization and polarimetric null columns can be set to zeros, for example: 

```
# Text file has columns wavelength, intensity, uncertainty
emptyarr = np.zeros_like(wavelength) # Arrays need to have the same length
spf_spectrum = pol.Spectrum(wavelength,intensity,emptyarr,emptyarr,emptyarr,uncertainty)
```

Additional details can be found in our tutorials. The [How to use Spectrum objects](../Tutorials/10-SpectrumClass_Tutorial.ipynb) or [How to use the LSD objects](../Tutorials/11-HowToLSD_Tutorial.ipynb) tutorials may be good places to start.


## Available Converters

SpecpolFlow provides converter functions that automatically reformat the files from some common spectropolarimeters into the `.s` LibreESPRIT format.

The current version of SpecpolFlow includes converters for the following spectropolarimetric data files:

| Instrument | Data Reduction Pipeline | Function |
|:---:|:---:|:---:| 
| ESPaDOnS | UPENA | converters.espadons() |
| SPIRou   | APERO | converters.spirou()   |

Narval observations are already in the correct format for SpecpolFlow. ESPaDOnS observations acquired from the PolarBase archive are also in the correct format.

For more details on converting observations, please refer to the [tutorial on converting to .s files](../Tutorials/1-ConvertToSFiles_Tutorial.ipynb). 

## Using SpecpolFlow in Spectroscopic Analyses 

SpecpolFlow can be used for any kind of spectroscopic data that fits the format described above. It has functions for 

* Normalizing spectra
* Coadding observations
* Merging spectral orders
* Isolating wavelength ranges or individual spectral lines
* Calculating bisectors, radial velocities, and equivalent widths

For more details, please refer to our suite of tutorials, especially the [Spectrum Class tutorial](../Tutorials/10-SpectrumClass_Tutorial.ipynb)
and [LSD Class tutorial](../Tutorials/5-LSDClass_Tutorial.ipynb).

If you have additional ideas about how to use SpecpolFlow or modify its existing tools, please see our [guidelines for contributing to SpecpolFlow](../Contributing/Contributing.md).
