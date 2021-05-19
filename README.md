[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# Te_HF_Conversion

Conversion of the elastic thickness of the lithosphere to the planetary heat flow..

## Description

**Te_HF_Conversion** is a simple code that allows to convert the elastic thickness of the lithosphere to planetary heat flow given several input parameters, including crustal thickness, strain rate, or radiogenic heating. The model makes use of the equating bending moment approach of [McNutt (1984)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB089iB13p11180) and has been used in [Broquet et al (2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL086746).

### Benchmarks
Heat flow calculations have been benchmarker to various published numbers in the litterature (e.g., [Solomon & Head, 1990](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB095iB07p11073). 

## Methods
`..` 

## Example scripts
`..` 

## How to install and run Displacement_strain_planet
If you would like to modify the source code, download the Displacement_strain_planet repository and install using pip (or pip3 depending on your installation).
```bash
    git clone https://github.com/AB-Ares/Te_HF_Conversion.git
    cd Te_HF_Conversion/
    pip install .
```
Alternatively, you can install Displacement-strain-planet via pip
```bash
   pip install Te-HF-Conversion
```

## To run the example scripts
```bash
    cd examples
    jupyter notebook Run_demo.ipynb
    python Mars_crust_displacement.py 
```

## Author
[Adrien Broquet](https://www.oca.eu/fr/adrien-broquet) (adrien.broquet@oca.eu)
