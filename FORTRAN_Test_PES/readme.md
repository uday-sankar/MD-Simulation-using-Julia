# User Notes
# Usage
To use this FORTRAN PES use the appropriate compiled fortran library.
Download this folder change the path variables for the PES in ```H2CO_triplet.jl```.
Use:
- ```CH2O_PIPNN_GNU_x86.s0```: Linux system
- ```CH2O_PIPNN.dylib```: for mac

Change ```LIBCH2O ``` and ```PES_PATH``` in H2CO_triplet.jl with the appropriate variable names in your system

You can also complie the FORTRAN file ```CH2O_PIPNN.f90``` as a shared library using gfortran for your OS.
## Files
The folder CH2O contain the weights and biases for the PIPNN.
This fodler is needed for using the PES

## Reference

This FORTRAN PES is not my own.
The PES was obtained from the chempotpy library which used a PIP-NN from another publication:
PIP NN: https://pubs.rsc.org/en/content/articlelanding/2017/cp/c7cp04578f, 
Li, Jun, Changjian Xie, and Hua Guo. "Kinetics and dynamics of the C (3 P)+ H 2 O reaction on a full-dimensional accurate triplet state potential energy surface." Physical Chemistry Chemical Physics 19, no. 34 (2017): 23280-23288.

Chempotpy: https://pubs.acs.org/doi/full/10.1021/acs.jpca.3c05Shu, 
Yinan, Zoltan Varga, Dayou Zhang, and Donald G. Truhlar. "ChemPotPy: A Python Library for Analytic Representations of Potential Energy Surfaces and Diabatic Potential Energy Matrices." The Journal of Physical Chemistry A 127, no. 45 (2023): 9635-9640.899,