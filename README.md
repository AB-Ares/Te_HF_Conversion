[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/368906695.svg)](https://zenodo.org/badge/latestdoi/368906695)

# Te_HF_Conversion

Elastic thickness of the lithosphere, yield strength envelope, and heat flow calculations.

## Description

**Te_HF_Conversion (TeHF)** is a simple code that allows to:

(1) Convert the elastic thickness of the lithosphere to planetary heat flow (and a yield strength envelope) given several input parameters including crustal thickness, strain rate, or radiogenic heating.

(2) Retrieve the elastic thickness (or mechanical thickness, and a yield strength envelope) of the lithosphere based on a temperature profile. 

The model makes use of the equating bending moment approach of [McNutt (1984)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB089iB13p11180) and has been used in [Broquet et al (2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL086746),[2021)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020JE006730), or [Maia & Wieczorek (2022)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021JE007004).

### Benchmarks
Heat flow calculations have been benchmarked to various studies (e.g., [McNutt, 1984](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB089iB13p11180) or [Solomon & Head, 1990](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB095iB07p11073)). I thank [Julia Maia](https://www.oca.eu/fr/julia-maia) for having performed some further benchmarks to the literature.

## Methods
`Curv_Moment`  Determine the bending moment given the input yield
strength envelope and curvature.

`Conversion_Te_HF`  Determine the surface, crustal, and mantle heat flows, mechanical thickness, and thermal gradients from input rheology and elastic parameters.

`Conversion_Tprofile_Te`  Determine yield strength envelope, mechanical thickness given the input temperature profile. Elastic thickness will also be output from the assumed plate curvature.

`Brittle_Strength` Compute the brittle part of the yield strength envelope based on planetary constants.

## Example scripts
`Mars_YSE`  Determine the surface, crust, and mantle heat flows for a given elastic thickness on Mars assuming a wet rheology for the diabase crust and olivine mantle. Plot the yield strength envelope.

`Mars_T_profile_Te`  Determine the yield strength envelope, mechanical and elastic thickness for an assumed temperature profile that is linear but has a cosine tapered temperature anomaly between 20 and 50 km depth.

`Venus_YSE`  Determine the surface, crust, and mantle heat flows for a given elastic thickness on Venus assuming a dry rheology for the diabase crust and olivine mantle. Plot the yield strength envelope.

## How to install and run Te_HF_Conversion
If you would like to modify the source code, download the Displacement_strain_planet repository and install using pip (or pip3 depending on your installation).
```bash
    git clone https://github.com/AB-Ares/Te_HF_Conversion.git
    cd Te_HF_Conversion/
    pip install .
```

## To run the example scripts
```bash
    cd examples
    python Mars_YSE.py 
    python Venus_YSE.py 
```

## Author
[Adrien Broquet](https://ab-ares.github.io/website/) (adrien.broquet@dlr.de)

## Cite
If you use the package, please cite the latest release as:
Adrien Broquet. AB-Ares/Te_HF_Conversion: 0.2.3 (Version 0.2.3). Zenodo. http://doi.org/10.5281/zenodo.4973893

This package was created for [Broquet et al., 2020](https://doi.org/10.1029/2019GL086746), which you can also cite as
```
    @article{Broquet2020,
    author = {Broquet, A. and Wieczorek, M. A. and Fa, W.},
    title = {Flexure of the Lithosphere Beneath the North Polar Cap of Mars: Implications for Ice Composition and Heat Flow},
    journal = {Geophysical Research Letters},
    volume = {47},
    number = {5},
    doi = {https://doi.org/10.1029/2019GL086746},
    year = {2020}}
```
