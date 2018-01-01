# thin_film_memetic_optimization

<img src=https://user-images.githubusercontent.com/34921690/34465238-55a1b1fc-ee5a-11e7-8b5c-13eb5e8a762a.png width="600" height="300" />

# Contents
This package can be used to design broadband thin-film spectral filters given a list of dispersive materials

* **algorithm:** This is a folder that contains the main mechanics of the memetic algorithm and its helper functions

* **dispersion data:** This is a folder that contains the dispersive refractive index of each material. The refractive indices of each material is stored in .xlsx format, where the first column specifies the wavelength (in microns), the second (third) column specifies the real (imaginary) part of the refractive indices

* **disp_generate_example.m:** This is an example script that shows the syntax of specifying the wavelengths of interest and interpolating the material's refractive index properties

* **optimization_example_cooler.m:** This is an example script that optimizes the reflection/emissivity spectrum of a radiative cooling device

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sidenote: [This article] discusses how a radiative cooling device works

# How to run the example optimization
1. Run "disp_generate_example.m" **once** to generate the refractive indices of materials in the wavelengths of interest
2. Run "optimization_example_cooler.m" to see the optimization of the spectrum of a radiative cooling device

# How to run your own optimizations
1. (Optional) Place your own materials' refractive indices in ".xlsx" format inside the "dispersion data" folder. First column is wavelength in microns, and second (third) column is the real (imaginary) part of the refractive index.
2. Edit "disp_generate_example.m" to specify the wavelengths that you are interested in and the materials you are interested in using. 
3. Edit "optimization_example_cooler" to specify:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Your target reflection spectrum*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Incidence and substrate materials*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *List of materials that you are using in your multi-layer thin film*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Fixed material layers between substrate and multi-layer thin film*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Angle of incidence*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Number of thin film layers*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Memetic algorithm parameters*  

# How to cite this work
The details of this algorithm is published [here in ACS Photonics]  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Article DOI: 10.1021/acsphotonics.7b01136  

[This article]: https://www.nature.com/articles/nature13883
[here in ACS Photonics]: http://pubs.acs.org/doi/abs/10.1021/acsphotonics.7b01136
 
