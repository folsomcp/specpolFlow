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
Spectropolarimetry is a primary investigative method in stellar astrophysics for interpreting the characteristics and geometry of aspherical systems. It is particularly useful for tracing the influence of magnetic fields on a host star and its environment. Maintaining modern, flexible, and accessible computational tools that enable spectropolarimetric studies is thus essential. The `SpecpolFlow` package is a new, completely Pythonic workflow for stellar spectropolarimetry. Its suite of tools provides a user-friendly interface for analyzing stellar spectra from an assortment of instruments and telescopes. `SpecpolFlow` contains tools for spectral normalization and visualization, the extraction of Least-Squares Deconvolution (LSD) profiles, the generation and optimization of line masks for LSD analyses, and the calculation of longitudinal magnetic field measurements from the LSD profiles. The `SpecpolFlow` website also includes an array of tutorials that guide users through standard analytic practices and demonstrate the capability of the software.`SpecpolFlow` is distributed as a free, open-source package that is licensed under the GNU General Public License v. 2.0, with fully documented tools (via an API) which are actively maintained by a dedicated team of contributors.

# Statement of need

Spectropolarimetry is an essential analytic technique in stellar astrophysics that is used to study the properties and geometry of a source and its environment. In the presence of symmetry-breaking phenomena such as stellar magnetic fields, spectropolarimetry can illuminate fundamental characteristics of the observed system. Magnetic fields are present in non-degenerate stars across the HR diagram, at every spectral classification and evolutionary stage [@DonatiLandstreet2009; @Mestel2012]. Cool stars (spectral types F, G, K, and M) have a theoretical magnetic incidence of 100%, and display a wide range of observed field strengths, driven by variations in their internal dynamos [@Reiners2012]. Approximately 10% of hot stars (spectral types O, B, and A), which make up less than 1% of main sequence stars, harbour strong magnetic fields [@Grunhut2017; @Sikora2019]. Off the main sequence, magnetic fields are found around evolved stars and compact stellar remnants such as white dwarfs and pulsars [@Ferrario2015]. Spectropolarimetry is a powerful tool for characterizing the strength, orientation, and topology of these magnetic fields; thus, the maintenance and dissemination of computational tools enabling spectropolarimetry is a key element of progressive research.

Spectropolarimetric studies of stellar magnetic fields typically leverage the splitting of individual spectral lines due to the Zeeman Effect. The Zeeman split components of the line are polarized, and they are shifted in wavelength proportional to the field strength. The longitudinal magnetic field ($\langle B_z \rangle$), the component of the vector field along the line-of-sight, will produce circularly polarized light [e.g., @Landstreet2015]. These shifts in wavelength are typically undetectable in the total intensity spectra, except for very strong fields; however the changes in polarization across a line provide a much more sensitive, unambiguous diagnostic. In practical observations, the polarization signal of the Zeeman effect in an individual line is often lost in the noise, thus it is important to efficiently combine information from many spectral lines. The Least-Squares Deconvolution [LSD, @Donati1997; @Kochukhov2010] approach is the most widely used method for detecting these polarization signatures. LSD is a multi-line technique, similar to cross-correlation, which produces a pseudo-average line profile at increased signal-to-noise relative to a single line. Spectral lines are selected using a "line mask" tailored to the star. Lines in the observation are modeled assuming they have a common shape and differ by individual factors called the "weight," which is essentially the convolution of a pseudo-average with a set of delta functions of different wavelength and amplitude. This model is inverted using a linear least-squares procedure, accounting for the uncertainties on observed spectral pixels, to derive the common line shape (the "LSD profile"). The estimate of $\langle B_z \rangle$ can then be computed from the resulting circular polarization and intensity profiles. The modeling of the variation of $\langle B_z \rangle$ as the star rotates enables further characterizations, such as determining the stellar rotational period and simple models of the magnetic topology. 

Several successful programs supporting spectropolarimetric analyses exist in the literature, although they are often not open source and are poorly documented. The original LSD code by @Donati1997 is efficient and effective, written in C, but it is both proprietary and undocumented.  A more recent program, `iLSD,` is an Interactive Data Language (IDL&reg;) interface around a Fortran core implementing LSD [@Kochukhov2010]. This code includes additional features such as the reconstruction of multiple line profiles from the same spectrum. It also includes a first-order Tikhonov regularization function as a way to minimize the degradation of the profiles' quality due to the impact of numerical noise in the deconvolution. However `iLSD` is also proprietary and has very limited documentation. Some additional support codes (e.g., to read and write LSD profiles, to create line masks, to combine observations) have been written in a variety of languages and passed down from person to person. While these codes are scientifically relevant, they also can present a series of meaningful challenges to prospective users. Legacy codes often lack thorough documentation (even the code itself can be opaque) and they typically do not come with a detailed users manual or set of tutorials. Those that lack a dedicated package distribution or active maintenance may have depreciated functionality, non-standard implementations, or other "quirks" that result from minimal version control or a limited distribution network. Software developed in languages like IDL&reg; may require payment of a subscription fee to obtain the necessary licensed permissions before the software can even be used. Factors like these negatively impact accessibility, particularly for students and early-career researchers, and can serve as a significant barrier to learning or to research productivity.

# Overview of SpecpolFlow

The `SpecpolFlow` package is a modernized, Pythonic revitalization of the computational tools that have preceded it, unifying the key elements of legacy codes into a single, accessible workflow. The software is open source, well documented, and the code itself is extensively commented and designed to be readable. `SpecpolFlow` produces consistent results with previous proprietary codes implementing similar algorithms. It produces consistent LSD profiles with the code of @Donati1997 and `iLSD`, and it produces consistent $\langle B_z \rangle$ values with the code of @Wade2000.

`SpecpolFlow` provides a unified toolset, with an ensemble of key Python classes and functions that:

- Convert observed spectra into a common file format
- Continuum normalize spectra, with an interactive graphical interface
- Calculate spectral activity indices
- Generate line masks for LSD from lists of atomic transitions
- Clean line masks to remove problem lines or regions
- Calculate LSD profiles
- Calculate radial velocity
- Calculate $\langle B_z \rangle$ values

Because the specpolFlow package has a pythonic, object-oriented design, it is flexible and versatile, such that users can build their own custom workflow from the available specpolFlow classes and functions. This can also significantly enhance scientific archiving and reproducability when combined with python notebooks and database packages such as [pandas](https://pandas.pydata.org/). This is especially empowering for researchers aiming to train and involve students in short term projects. 

The continuum normalization routine follows the algorithm briefly described in @Folsom2008 and @Folsom2013, and implements a graphical interface and plotting with the Tkinter and matplotlib packages. The activity indices are calculated using the method and calibration coefficients of @Marsden2014. Line masks can be generated from line lists in the Vienna Atomic Line Database [@Ryabchikova2015] ``extract stellar long'' format, and if necessary, effective Landé factors can be estimated in LS, J$_1$J$_2$, and J$_1$K coupling schemes [@Martin1978; @LandiAndLandolfi2004]. The LSD calculation follows the method of @Donati1997, and more specifically the details from @Kochukhov2010. It relies on [numpy](https://numpy.org/) and makes careful use of numpy's sparse arrays for efficiency. Radial velocities are calculated by fitting a Gaussian to an LSD profile by default, although calculating first moments is also supported. The $\langle B_z \rangle$ calculation uses the first moment technique applied to LSD profiles [@Rees1979; @Donati1997; @Kochukhov2010].

`SpecpolFlow`'s [website](https://folsomcp.github.io/specpolFlow/) includes an extensive set of tutorials that explain the use of the package and demonstrate some common analytic cases. The tutorials also include suggestions for using the workflow in a collaborative environment (using [Google Colab](https://colab.research.google.com/)) and automation methods for large datasets. These tutorials, in the form of [Jupyter Notebooks](https://jupyter.org/), can even be used within a classroom or workshop setting.

Due to its flexibility and versatility, `SpecpolFlow` has already been used in several projects. The python package [`pyRaven`](https://veropetit.github.io/pyRaven/) incorporates `SpecpolFlow` to determine the magnetic dipolar field upper limits, using the Bayesian method from @Petit2012. The package [`ZDIpy`](https://github.com/folsomcp/ZDIpy) [@Folsom2018] uses results from `SpecpolFlow` to construct magnetic maps using Zeeman-Doppler Imaging. It has also been used for a new spectropolarimetric analysis of the most recently discovered magnetic O-type star HD 54879 [@Erba2024] and for detecting an exoplanet candidate around a young K-type star [@Zaire2024]. 

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


