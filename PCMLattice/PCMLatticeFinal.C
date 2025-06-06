/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Reference
    Dominik Wagner (dominik.wagner1@rwth-aachen.de)
    Aachen Institute of Structural Mechanics and Lightweight Design (SLA)

Application
    PCMLatticeFinal

Description
    Transient solver for incompressible convection diffusion phase change. The solver is based on
    buoyantBoussinesqPimpleFoam. Convection is induced by Boussinesq approximation.
    Integration of porous media lattice structure via permeability. Assumption of local thermal equilibrium using effective material properties.

    Solver is derived from CoMeTFoam Solver by Kai Schueller (schueller@aices.rwth-aachen.de)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicTransportModel.H"
//#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvModels.H"
#include "fvConstraints.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvModels.H"
    #include "createFvConstraints.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // read some control values
    const dictionary& phasechangeDict = mesh.solutionDict().subDict("PHASECHANGE");
    Foam::label minTiter(readLabel(phasechangeDict.lookup("minTiter")));
    Foam::label maxTiter(readLabel(phasechangeDict.lookup("maxTiter")));
    scalar gamma_L_Tol(readScalar(phasechangeDict.lookup("gamma_L_Tol")));
    scalar relaxGamma_L(readScalar(phasechangeDict.lookup("relaxGamma_L")));

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
