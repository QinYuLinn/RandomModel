import numpy as np
from TenFoldHamiltonian import *#import the class defined for the symmetry operator
dim=0#real space dimension
lattice_size=[1]#site numbers in each dimesion
spin=3/2#internal degree of freedom in each lattice site
symmetry_vec=[0,0,0]#the symmetry class vector,T,C,S,0 represent without this symmetry

SymOpe=TenFoldHamiltonian(dim,lattice_size,spin,symmetry_vec)#construct the symmetry operator
SymOpe.ConstructSymmetryOperators()
SymOpe.ConstructTenfoldHamiltonian()
SymOpe.CalculateTopologicalInvariant()