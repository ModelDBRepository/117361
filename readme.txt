Readme for the model code associated with the paper:

Stewart RD, Bair W. (2009)
Spiking neural network simulation: numerical integration with the Parker-Sochacki method.
J Comput Neurosci
ISSN:	0929-5313 (Print) 1573-6873 (Online)
DOI:	10.1007/s10827-008-0131-5
Springer Netherlands
SpringerLink Date: Saturday, January 17, 2009

The article is freely available through Springer's Open Access program. It can be downloaded from here:
http://www.springerlink.com/content/8q431074727r41h2/

The code was uploaded to the ModelDB database on 13/02/09. It is divided into three folders:

1) mex - contains all the c code, plus compiled mex files of the top-level simulation functions iz_ps.c 
and tm_ps.c (both 64-bit Windows and 64-bit Linux versions). Instructions for compiling are given in the 
comments at the top of the .c files.

2) iz - contains the matlab code for the Izhikevich model simulations.
iz_bench.m is the top level script for running all experiments.

3) tm - contains the matlab code for the Traub-Miles HH neuron simulations.
tm_bench.m is the top level script for running all experiments.

Also included in both iz and tm folders are .mat files with results from the injection current simulations. 
The Izhikevich recurrent network simulation results take up GBs of space and have not been uploaded. 
To generate them locally, you need to run iz_bench and select simulation type 4.