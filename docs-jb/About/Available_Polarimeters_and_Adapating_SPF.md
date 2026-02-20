# Compatible Spectropolarimeters and Spectrographs  

SpecpolFlow is extensible to various spectropolarimeters and even spectrographs. 

To use SpecpolFlow's main functions, your data's structure must follow the LibreESPRIT format.

## Available Converters
We can develop converters that reformat observations to provide a smooth workflow for users.  

The current version of SpecpolFlow includes converters for the following spectropolarimetric data files:

| Instrument | Data Reduction Pipeline | Function |
|:---:|:---:|:---:| 
| ESPaDOnS | UPENA | converters.espadons() |
| SPIRou   | APERO | converters.spirou()   |

Narval observations are already in the correct format for SpecpolFlow. ESPaDOnS observations acquired from the PolarBase archive are also in the correct format.
For more details on converting observations, please refer to the [tutorial on converting to .s files](../Tutorials/1-ConvertToSFiles_Tutorial.ipynb). 

## Using SpecpolFlow in Spectroscopic Analyses 

SpecpolFlow can be used for spectroscopic data and is not exclusively limited to stellar astronomy. 
As long as the data contain arrays of flux as a function of wavelength, SpecpolFlow can reasonably be extended to different kinds of spectroscopic datasets,
such as synthetic spectra, space-based spectroscopy, UV bands, etc.. Some of SpecpolFlows's functions can be used to: 

* Manage observational data
* Normalize spectra
* Coadd observations
* Isolate wavelength ranges or spectra lines 
* Calculate bisectors   
* Calculate radial velocities

For more details refer to our suite of tutorials, especially the [Spectrum Class tutorial](../Tutorials/10-SpectrumClass_Tutorial.ipynb)
and [LSD Class tutorial](../Tutorials/5-LSDClass_Tutorial.ipynb)