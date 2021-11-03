[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/368906695.svg)](https://zenodo.org/badge/latestdoi/368906695)

# Te_HF_Conversion

Conversion of the elastic thickness of the lithosphere to the planetary heat flow.

## Description

**Te_HF_Conversion** is a simple code that allows to convert the elastic thickness of the lithosphere to planetary heat flow given several input parameters, including crustal thickness, strain rate, or radiogenic heating. The model makes use of the equating bending moment approach of [McNutt (1984)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB089iB13p11180) and has been used in [Broquet et al (2020)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL086746).

### Benchmarks
Heat flow calculations have been benchmarked to various studies (e.g., [McNutt, 1984](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB089iB13p11180) or [Solomon & Head, 1990](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB095iB07p11073)). I thank [Julia Maia](https://www.oca.eu/fr/julia-maia) for having performed some further benchmarks to the literature.

### Caution
This code doesn't properly account for crust/mantle decoupling when estimating elastic strengths (see [Burov & Diament, 1984](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/94JB02770)). The decoupled mantle doesn't have its own elastic core. This will be fixed soon.

## Methods
`Curv_Moment`  Determine the bending moment given the input yield
strength envelope and curvature.

`CalcYSE_HF`  Determine the surface, crustal, and mantle heat flows, mechanical thickness, and thermal gradients from input rheology and elastic parameters.

## Example scripts
`Mars_YSE`  Determine the surface, crust, and mantle heat flows for a given elastic thickness on Mars assuming a wet rheology for the diabase crust and olivine mantle. Plot the yield strength envelope.

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
[Adrien Broquet](https://www.oca.eu/fr/adrien-broquet) (adrien.broquet@oca.eu)

## Cite
You can cite the latest release of the package as:
Adrien Broquet. AB-Ares/Te_HF_Conversion: 0.1.1 (Version 0.1.1). Zenodo. http://doi.org/10.5281/zenodo.4973893
