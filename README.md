This is a code to generate IC's with Eddingtons formula

How to compile?
===========

Simply just run "make". Make sure "gsl" (GNU Scientific Library) is installed/configured.
 

How to run the code
===========

To compile and run the code do the following:

	make
	./EddIso Output.bin 1000000

This will create a Gadget2-file (Output.bin) with 1000000 particles. Note that the block with potentials are omitted from the output file.


The parameters
===========


Some of the parameters that can be modified are listed here:

The density profile is determined by the following lines (in EddIso.c):

	//DENSITY PROFILE: (inner slope, outer slope, transition...)
	#define ALPHA (1.0)
	#define BETA (4.0)
	#define GAMMA (1.0)

(With the given parameters it is a Hernquist profile)

The density profile is

	Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - alpha) / gamma);

Here Rho0 is a constant, and a is the scalelength.


With the following lines (in EddIso.c) Rho0 and the scale_length and the virial radius can be set

	Cl->Rho0 = 1.0/2.0/PI;
	Cl->scale_length = 1.0;
	Cl->rvir = 50.;

The gravitational constant can be set with

	#define G (1.0) // Gravitational constant
