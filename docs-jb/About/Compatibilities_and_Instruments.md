# Compatibilities

SpecpolFlow is generally compatible with data in the [LibreESPRIT](https://www.cfht.hawaii.edu/Instruments/Spectroscopy/Espadons/Espadons_esprit.html) format. 
For spectropolarimetric data, this is usually a text file with columns for wavelength, intensity, polarization, two polarimetric nulls, and uncertainties.

For intensity spectra (with no polarization) this is a text tile with column for wavelength, intensity, and uncertainties.

Additional details can be found in the [converting to .s files](../Tutorials/1-ConvertToSFiles_Tutorial.ipynb) and [spectrum objects](../Tutorials/10-SpectrumClass_Tutorial.ipynb) tutorials, which use the [read_spectrum](read_spectrum) function.  


## Available Converters

SpecpolFlow provides converter functions that automatically reformat the files from some common spectropolarimeters into the `.s` LibreESPRIT format.

The current version of SpecpolFlow includes converters for the following spectropolarimetric data files:

| Instrument | Data Reduction Pipeline | Function |
|:---:|:---:|:---:| 
| ESPaDOnS | UPENA | converters.espadons() |
| SPIRou   | APERO | converters.spirou()   |

Narval observations are already in the correct format for SpecpolFlow. ESPaDOnS observations acquired from the PolarBase archive are also in the correct format.

For more details on converting observations, please refer to the [tutorial on converting to .s files](../Tutorials/1-ConvertToSFiles_Tutorial.ipynb).

If you are interested in developing a converter for another instrument, you are welcome to contact us, or check the [contributing](../Contributing/Contributing.md) page.


## Using SpecpolFlow in Spectroscopic Analyses 

SpecpolFlow can be used for any kind of spectroscopic data that fits the format described above. The SpecpolFlow tools will generally work with echelle spectra and single order spectra. Although observations with no formal uncertainty may cause problems in some functions, and result in bad uncertainty estimates. It has functions for: 

* Normalizing spectra
* Coadding observations
* Merging spectral orders
* Isolating wavelength ranges or individual spectral lines
* Calculating bisectors, radial velocities, and equivalent widths

For more details, please refer to our suite of tutorials, especially the [Spectrum Class tutorial](../Tutorials/10-SpectrumClass_Tutorial.ipynb)
and [LSD Class tutorial](../Tutorials/5-LSDClass_Tutorial.ipynb).

If you have additional ideas about how to use SpecpolFlow or modify its existing tools, please see our [guidelines for contributing to SpecpolFlow](../Contributing/Contributing.md).
