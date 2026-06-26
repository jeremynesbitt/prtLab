# prtLab 

prtLab is a Matlab-primary toolkit to do true Polarization Ray Tracing of simple object such as waveplates.  It not a commercial solution but is usable as a research tool.

# Polarization Ray Tracing

This project was inspired by the excellent book Polarized Light and Optical Systems by Chipman, Lam and Young (CRC Press).

This book has numerous examples with great plots to explain the framework.  I started out for fun to reproduce some of the plots (see https://github.com/jeremynesbitt/PolarizedLightAndOpticalSystems), but after getting deeper into the core of the book I did not see much published code that used the P matrix framework.  So I decided to separate it out.

# Project Structure

The best way to engage with the project is to look at the examples, and build off of them.  See the examples and the docs tab.

Matlab source code is in the prtLab folder.  There is a python version, which is mainly an AI generated conversion of the matlab code.  

# Installation

For matlab, just copy the prtLab folder to your local directory and add it to your path.  It mostly is okay with the base package; only the asphere ray tracing requires the optimization toolbox at present.

