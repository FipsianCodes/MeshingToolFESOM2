
A meshing tool to refine the PI mesh for FESOM2 simulations, provided by the Alfred-Wegener-Institute. The mesh is obtainable:
https://gitlab.awi.de/fesom/

The meshing tool is written in Python and is based on the packages:
  numpy
  dataclasses 
  matplotlib

An additional option is provided to write out gridfiles, that are needed to perform conservative interpolation methods. The PI mesh 
and the refinement FPI were used to apply a Micro-Macro Parareal algorithm to FESOM2, see:

https://arxiv.org/abs/2306.17269

WARNING:
The plotting option was implemented to test the refinement procedure for very small sub sets of the PI domain. In case of a FPI mesh
vizualisation the performance significantly slows down.
