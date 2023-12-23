# PCLatFoam
OpenFoam solver for the simulation of the thermal and fluid flow behaviour of Phase Change Material / Lattice structure composites. The solver is an extension and readaptation the OpenFOAMv9 of the CoMeTFoam solver (https://github.com/geo-fluid-dynamics/CoMeTFoam).
The solver includes a Darcy-Brinkman term to account for the fluid flow within the porous medium (lattice structure). Local Thermal Equilibrium (LTE) is used for the energy equation. However, the convective terms see corrected fluid properties with respect to the effective ones.
