---
title: 'SpecpolFlow: a new software package for spectropolarimetry using Python'
tags:
  - Spectropolarimetry
  - Stellar magnetic fields
  - Astronomy software
  - Open source software
  - Astronomy data visualization
authors:
  - name: Colin P. Folsom^[Co-first author]
    orcid: 0000-0002-9023-7890
    corresponding: true
    affiliation: 1
  - name: Christiana Erba^[Co-first author]
    orcid: 0000-0003-1299-8878
    affiliation: 2
  - name: Veronique Petit
    orcid: 0000-0002-5633-7548
    affiliation: 3
  - name: Shaquann Seadrow
    orcid: 0009-0002-0308-2497
    affiliation: 3
  - name: Patrick Stanley
    orcid: 0000-0002-0378-0140
    affiliation: 3 
  - name: Tali Natan
    orcid: 0000-0002-7703-6701
    affiliation: 3
  - name: Bonnie Zaire
    orcid: 0000-0002-9328-9530
    affiliation: 4
  - name: Mary E. Oksala
    orcid: 0000-0003-2580-1464
    affiliation: "5, 6"
  - name: Federico Villadiego-Forero
    #orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Robin Moore
    #orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Marisol Catalan Olais
    orcid: 0009-0006-2442-6235
    affiliation: 3    
   
affiliations:
 - name: Tartu Observatory, University of Tartu, Observatooriumi 1, 61602, Toravere, Estonia
   index: 1
 - name: Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
   index: 2
 - name: Department of Physics and Astronomy, Bartol Research Institute, University of Delaware, 19716, Newark, DE, USA
   index: 3
 - name: Universidade Federal de Minas Gerais, Belo Horizonte, MG, 31270-901, Brazil
   index: 4
 - name: Department of Physics, California Lutheran University, 60 West Olsen Road, 91360, Thousand Oaks, CA, USA
   index: 5
 - name: LESIA, Observatoire de Paris, PSL University, CNRS, Sorbonne Université, Université Paris Cité, 5 place Jules Janssen, 92195 Meudon, France
   index: 6
 
date: 25 October 2024
bibliography: paper.bib

---

# Summary

Spectropolarimetry, the observation of polarization and intensity as a function of wavelength, is a powerful tool in stellar astrophysics. It is particularly useful for characterizing the distribution of circumstellar material and for tracing the influence of magnetic fields on a host star and its environment. Maintaining modern, flexible, and accessible computational tools that enable spectropolarimetric studies is thus essential. The `SpecpolFlow` package is a new, completely Pythonic workflow for analyzing stellar spectropolarimetric observations. Its suite of tools provides a user-friendly interface for working with data from an assortment of instruments and telescopes. `SpecpolFlow` contains tools for spectral normalization and visualization, the extraction of Least-Squares Deconvolution (LSD) profiles, the generation and optimization of line masks for LSD analyses, and the calculation of longitudinal magnetic field measurements from the LSD profiles. It also provides Python classes for the manipulation of spectropolarimetric products. The `SpecpolFlow` website includes an array of tutorials that guide users through common analytic analysis using the software. `SpecpolFlow` is distributed as a free, open-source package, with fully documented tools (via an API and command line interface) which are actively maintained by a team of contributors. 

# Statement of need

Spectropolarimetry is an essential observational technique in stellar astrophysics, which is used to study the surface properties of stars and their environments. Symmetry-breaking phenomena like stellar magnetic fields leave an imprint on polarized spectra.  Magnetic fields are present in most classes of stars throughout their evolution [@DonatiLandstreet2009; @Mestel2012]. Cool stars (with convective envelopes) have a theoretical magnetic incidence of 100%, and display a wide range of observed field strengths, driven by variations in their internal dynamos [@Reiners2012]. Nearly 10% of hot stars (with radiative envelopes) also harbour strong magnetic fields [@Grunhut2017; @Sikora2019]. Furthermore, magnetic fields are found in evolved giants and compact stellar remnants such as white dwarfs and pulsars [@Ferrario2015]. Spectropolarimetry is a valuable tool for characterizing the strength, orientation, and topology of these fields. Therefore the maintenance and dissemination of computational tools enabling this technique are a key element of research.

Spectropolarimetric studies of stellar magnetic fields typically leverage the splitting of spectral lines due to the Zeeman effect [e.g., @Landstreet2015]. The Zeeman split components of a line are polarized and shifted in wavelength proportional to field strength. These shifts in wavelength are typically undetectable due to other line broadening processes (except for very strong fields); however the changes in polarization remain detectable. In practical observations, within an individual line, this polarization signal is often below the noise level, thus it is important to combine information from many spectral lines. The Least-Squares Deconvolution [LSD, @Donati1997; @Kochukhov2010] approach is the most widely used method for detecting such polarization signatures. LSD is a multi-line technique, similar to cross-correlation, which produces a pseudo-average line profile at increased signal-to-noise. LSD models the spectrum as the convolution of a set of delta functions (the 'line mask') with a common line shape (the 'LSD profile'). This model is fit to an observation, using a weighted linear least-squares approach, to derive the LSD profile. An estimate of the surface averaged line-of-sight component of the field (the 'longitudinal field', $\langle B_z \rangle$) can then be computed from the circular polarization and intensity LSD profiles. Modelling the variation of $\langle B_z \rangle$ as the star rotates enables further characterization, such as determining the stellar rotational period and simple models of the magnetic topology.

Several successful programs supporting spectropolarimetric analyses exist in the literature, although they are often not open source and poorly documented. There is currently no publicly available software that provides the full toolset needed to analyze reduced spectropolarimetric observations and produce magnetic field measurements. The original LSD code by @Donati1997 is efficient and effective, written in C, but it is both proprietary and undocumented.  A more recent program, `iLSD,` is an Interactive Data Language (IDL) interface around a Fortran core implementing LSD [@Kochukhov2010]. This code includes additional features such as the reconstruction of multiple line profiles from the same spectrum, or an optional Tikhonov regularization of the LSD profiles. However `iLSD` is also proprietary and has very limited documentation. Some additional support codes (e.g., to read and write LSD profiles, to create line masks, or to combine observations) have been written in a variety of languages and passed down from person to person. While these codes are scientifically relevant, they often lack documentation, version control, active maintenance, and generally they are neither publicly available nor open source. Software developed in proprietary languages like IDL may require a subscription fee before the software can even be used. Factors like these negatively impact accessibility, particularly for students and early-career researchers, and can serve as a significant barrier to both learning and research reproducibility.

# Overview of SpecpolFlow

The `SpecpolFlow` package is a modernized, unified, Pythonic revitalization of the computational tools that have preceded it. The software is open source, well documented, and the code itself is extensively commented and designed to be readable. `SpecpolFlow` produces consistent results with previous proprietary codes implementing similar algorithms. It produces consistent LSD profiles with the code of @Donati1997 and `iLSD`, and it produces consistent $\langle B_z \rangle$ values with the code of @Wade2000.

`SpecpolFlow` provides a toolkit with an ensemble of Python functions that:

- Convert observed spectra into a common file format
- Continuum normalize spectra, with an interactive graphical interface
- Generate line masks for LSD from lists of atomic transitions
- Clean line masks interactively to remove problem lines or adjust line depths
- Calculate LSD profiles
- Calculate line bisectors
- Calculate radial velocities
- Calculate $\langle B_z \rangle$ values

The continuum normalization routine follows the algorithm briefly described in @Folsom2008 and @Folsom2013, and implements a graphical interface and plotting with the Tkinter and matplotlib packages. Line masks can be generated from line lists in the [Vienna Atomic Line Database](https://vald.astro.uu.se/) [@Ryabchikova2015] "extract stellar" "long" format. If necessary, effective Landé factors are estimated in LS, J$_1$J$_2$, and J$_1$K coupling schemes [@Martin1978; @LandiAndLandolfi2004]. The LSD calculation follows the method of @Donati1997, with details from @Kochukhov2010. It relies on [numpy](https://numpy.org/) and makes careful use of numpy's sparse arrays for efficiency.
The interactive line mask cleaning tool uses Tkinter and matplotlib for the interface. This tool can remove lines from the mask, and automatically fit line depths using a reference observed spectrum, following @Grunhut2017 with some optimizations. The line depth fitting routine inverts the linear least squares problem in LSD, and instead fits line depths given an observation and LSD profile. That LSD profile must be approximately correct, calculated using a mask (or part of a mask) with dominantly acceptable lines. Line depth fitting remains a weighted linear least squares problem, similar to the LSD calculation, and can be solved efficiently using spare matrix operations. Radial velocities are calculated by fitting a Gaussian to an LSD profile by default, although calculating from first moments is also supported. The $\langle B_z \rangle$ calculation uses the first moment technique applied to LSD profiles [@Rees1979; @Donati1997; @Kochukhov2010].

SpecpolFlow enables users to build their own custom workflow from the available classes and functions. This can significantly enhance scientific archiving and reproducibility when combined with Jupyter notebooks or database packages such as [pandas](https://pandas.pydata.org/). This is especially valuable for training or involving students in projects. To aid this, `SpecpolFlow`'s [website](https://folsomcp.github.io/specpolFlow/) includes an extensive set of tutorials that explain the use of the package and demonstrate some common analysis cases. The tutorials also include suggestions for using a collaborative environment ([Google Colab](https://colab.research.google.com/)) and automation methods for large datasets. These tutorials, in the form of [Jupyter Notebooks](https://jupyter.org/), can flexibly be used within a classroom or workshop setting.

`SpecpolFlow` has already been used in several projects. The python package [`pyRaven`](https://veropetit.github.io/pyRaven/) incorporates `SpecpolFlow` to determine the magnetic dipolar field upper limits, using the Bayesian method from @Petit2012. The package [`ZDIpy`](https://github.com/folsomcp/ZDIpy) [@Folsom2018] uses results from `SpecpolFlow` to construct magnetic maps using Zeeman-Doppler Imaging. It has also been used for a new spectropolarimetric analysis of the most recently discovered magnetic O-type star HD 54879 [@Erba2024] and for detecting an exoplanet candidate around a young K-type star [@Zaire2024]. 

# Acknowledgements
The authors gratefully acknowledge past contributions to the early development and testing of `SpecpolFlow`'s tools, including Dax Moraes, David Meleney Jr., and Gregg Wade.

This research was supported by the Munich Institute for Astro-, Particle and BioPhysics (MIAPbP), which is funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany´s Excellence Strategy – EXC-2094 – 390783311.

CPF gratefully acknowledges funding from the European Union's Horizon Europe research and innovation programme under grant agreement No. 101079231 (EXOHOST), and from the United Kingdom Research and Innovation (UKRI) Horizon Europe Guarantee Scheme (grant number 10051045). 

VP, SS, PS, and TN gratefully acknowledge support for this work from the National Science Foundation under Grant No. AST-2108455. 

MEO gratefully acknowledges support for this work from the National Science Foundation under Grant No. AST-2107871. 

SS gratefully acknowledges support for this work from the Delaware Space Grant College and Fellowship Program (NASA Grant 80NSSC20M0045).

BZ acknowledges funding from the CAPES-PrInt program (#88887.683070/2022-00 and #88887.802913/2023-00).

The `SpecpolFlow` team also thanks Ms. Tali Natan for her creative design of the official SpecpolFlow logo. 

# References

