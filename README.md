# prtLab 

prtLab is a Matlab-primary toolkit to do true Polarization Ray Tracing of simple object such as waveplates.  It not a commercial solution but is usable as a research tool.

# Polarization Ray Tracing

This project was inspired by the excellent book Polarized Light and Optical Systems by Chipman, Lam and Young (CRC Press, ISBN 9781498700566).

In working through some of the examples of the book, I did not find any prior open source work that used this framework.  So for fun, I decided to 

# Project Structure

The best way to engage with the project is to look at the examples, and build off of them.  See the examples and the docs folders.  The docs folder are .md files based on the .ipynb files which go through many of the function calls in a bit of detail.

Matlab source code is in the prtLab folder.  There is a python version, which is mainly an AI generated conversion of the matlab code.  

# Installation

For matlab, just copy the prtLab folder to your local directory and add it to your path.  It mostly is okay with the base package; only the asphere ray tracing requires the optimization toolbox at present.  

# Future Work

I would like to add some more examples, such as more crystal polarizers.  I also plan to add support for biaxial materials, since I punted on this initially.  Whether any more improvements happen after this I am not sure, but suggestions are welcome.

