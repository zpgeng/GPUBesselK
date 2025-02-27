# GPU-Acceleration of Modified Bessel Function of the Second Kind

This is the repo for the manuscript *GPU-Accelerated Modified Bessel Function of the Second Kind for Gaussian Processes* (https://arxiv.org/pdf/2502.00356)

File hierarchy:

- `eval_bessel.c`  The main program that is to be compiled and executed, which generates an array of BesselK function evaluations.
- `src` This folder contains the CUDA program of executing BesselK function.
- `obj` This folder contains object files.
- include This folder contains the header file that is to be used by eval_bessel.c.
