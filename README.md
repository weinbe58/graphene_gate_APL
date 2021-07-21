# graphene_gate_APL
This code is used in the calculations of

# Exact diagonalization (toy_model_ED.py)

In order to run the python code one needs to install [QuSpin](https://github.com/weinbe58/QuSpin) (installation instructions found on the QuSpin website) and should be executed with `Kondo_basis.py` in the same directory as `toy_model_ED.py` 

# DMRG calculation
The source code is all in this repo all that is required is lapack and blas installed in the standard location. The Makefile works for OS-X and linux. with clang and g++.



# extra file

## mappings/single_site_graphene.txt 
This file is used in both the ED and DMRG calculations. This file contains the mapping for the graphene flakes. 
