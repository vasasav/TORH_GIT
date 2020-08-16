# Toroidal Harmonics

This package provides python tools and bindings to work with toroidal harmonics

## Introduction

TORH is a part of the larger research programme aimed at creating software to study electromagnetic non-radiating configurations and Flying Doughnut pulses (Nature Mater. 15, 263 (2016) or https://rdcu.be/6I2l). In particular, one of the aims of the programme is to identify geometries (e.g. of dielectric and plasmonic particles) which could support strong non-radiating configurations. Current understanding of non-radiating configurations suggests that geometries that support them could often have toroidal shape. In search for these geometries it is therefore beneficial to work in toroidal coordinates (https://en.wikipedia.org/wiki/Toroidal_coordinates).

Unfortunately the numerical eco-system for working with toroidal coordinates is quite underdeveloped. The aim of TORH is therefore to provide tools, mostly in Python, for functional analysis in toroidal coordinates. This includes the effective numerical implementation of decompositions of scalar and vector fields into the so-called TORoidal Harmonics (TORH).

## Large-Scale Structure

The library is still under development. The current structure is:

* `DTORH` and `DTORH64` a 32- and 64-bit wrappers around the Fortran routines for computing associated Legendre polynomials of fractional order for x>1. Original code is by Segura and Gil (Comput. Phys. Comm. 124, 104 (2000)). The code is compiled into Dynamic Linked Libraries (`wrapDTORH.dll` and `wrapDTORH64.dll`) with wrappers in C. See `DTORH64/wrapDTORH.h` for function definitions. Due to packaging of the necessary fortran routines, currently the package runs on Windows only.

* `ToroidalHarmonics/DTORH.py` is provides Python (NumPy) wrapper around `wrapDTORH64.dll`, including handling of numerical overflows via standard Python exceptions (rather than Fortran error codes). 
    
* `ToroidalHarmonics/TorHarmRep.py` provides toolset for the decomposition of scalar functions, defined in 3D space, into toroidal harmonics (solutions of Laplace equation in toroidal coordinates). See `ToroidalHarmonics/DemoTorHarmRep.py` for demonstration.
    
* `ToroidalHarmonics/TorHarmVecRep.py` provides toolset for the decomposition of vector-valued functions, defined in 3D space, into vectorial toroidal harmonics. See `ToroidalHarmonics/testToroHarmVecRep.py` for demonstration. The theory for vectorial toroidal harmonics had to be developed as well, and is available in a summary form in `ToroidalHarmonics/ToroidalHarmonicsDefinition_Convention/ Toroidal_Harmonics_Def2.pdf`. 
