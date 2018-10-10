# RSBEM
regularized Stokeslet Boundary Element Method in MATLAB

# README #

### What is this repository for? ###

This repository contains a MATLAB implementation of a regularized Stokeslet Boundary Element Method (RSBEM).  It is heavily inspired by (1) but utilizes 2nd order curved triangular surface meshes and some performance improvements.  So far, it has primarily been used to model swimming monoflagellated bacteria (2) as well as swimming dinoflagellates (3) but effort has been made to keep the code general where possible, so other geometries and low-Reynolds number problems can be modeled in the future.

The Salome Platform (https://www.salome-platform.org/) has been used to generate geometries and meshes; python scripts related to this for curved rod bacteria and helical bacterial flagella are also included.

**This code is a work in progress.  It is not yet documented, easily useable, etc, but I am working toward this goal.**

(1) D. J. Smith, A boundary element regularised Stokeslet method applied to cilia and flagella-driven flow. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences. 465, 3605–3626 (2009).
(2) R. Schuech, T. Hoehfurtner, D. Smith, and S. Humphries.  Motile curved bacteria are Pareto-optimal.  in prep. for submission to Science
(3) L. T. Nielsen, R. Schuech, S. Humphries, D. Smith, and T. Kiørboe.  Hydrodynamics shed light on dinoflagellate ecology and evolution.  in prep. for submission to eLife


### How do I get set up? ###

All code (with notable exception of Python scripts for geometry / meshing in Salome) is in MATLAB.  The Statistics and Machine Learning Toolbox is needed as well as the Parallel Computing Toolbox for a number of parallelized routines (e.g. using parfor).  In addition, several computationally intensive functions (i.e. calls to "...mexed()" functions) can/should be compiled to C using the Coder Toolbox - this yields a factor of ~30 speed increase in prior tests.




### Contribution guidelines ###

Please send any comments, inquiries, etc to rudi.schuech@gmail.com.
