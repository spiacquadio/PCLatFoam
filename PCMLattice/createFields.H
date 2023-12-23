Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

#include "readTransportProperties.H"

//Calculating phase-fraction gamma_L
Info<< "Calculating volume fraction gamma_L\n" << endl;
volScalarField gamma_L
(
    IOobject
    (
        "gamma_L",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    max(min((T-T_S)/(T_L-T_S),scalar(1)),scalar(0))
);
gamma_L.write();

volScalarField Tcons("Tcons", (T_L-T_S)*gamma_L + T_S);


volTensorField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh
    ),
    gamma_L*(k_L) + (scalar(1) - gamma_L)*(k_S)
);


volScalarField cpeff
(
    IOobject
    (
        "cpeff",
        runTime.timeName(),
        mesh
    ),
    gamma_L*(cpeff_L)
            + (scalar(1) - gamma_L)*(cpeff_S)
);
cpeff.oldTime();

volScalarField cpPCM
(
    IOobject
    (
        "cpPCM",
        runTime.timeName(),
        mesh
    ),
    gamma_L*(cpPCM_L)
            + (scalar(1) - gamma_L)*(cpPCM_S)
);
cpPCM.oldTime();

volScalarField rhoeffmix
(
    IOobject
    (
        "rhoeffmix",
        runTime.timeName(),
        mesh
    ),
    gamma_L*(rhoeff_L)
    +(scalar(1)-gamma_L)*(rhoeff_S)
);

volScalarField rhoPCMmix
(
    IOobject
    (
        "rhoPCMmix",
        runTime.timeName(),
        mesh
    ),
    gamma_L*(rhoPCM_L)
    +(scalar(1)-gamma_L)*(rhoPCM_S)
);

volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh
    ),
    nua+(scalar(1)-gamma_L)*nun //for velocity 0
);

volScalarField nu_eff
(
    IOobject
    (
        "nu_eff",
        runTime.timeName(),
        mesh
    ),
    (nu_effL)+(scalar(1) - gamma_L)*(nu_effS)
);


volScalarField rhoeffmixCpeff
(
    IOobject
    (
        "rhoeffmix*cpeff",
        runTime.timeName(),
        mesh
    ),
    rhoeffmix*cpeff
);

surfaceScalarField rhoeffmixPhiCpeff
(
    IOobject
    (
        "rhoeffmix*phi*cpeff",
        runTime.timeName(),
        mesh
    ),
    phi*fvc::interpolate(cpeff)*fvc::interpolate(rhoeffmix)
);

surfaceScalarField rhoPCMmixPhiCpPCM
(
    IOobject
    (
        "rhoPCMmix*phi*cpPCM",
        runTime.timeName(),
        mesh
    ),
    phi*fvc::interpolate(cpPCM)*fvc::interpolate(rhoPCMmix)
);

surfaceScalarField rhoeffmixPhi
(
    IOobject
    (
        "rhoeffmix*phi",
        runTime.timeName(),
        mesh
    ),
    phi*fvc::interpolate(rhoeffmix)
);

surfaceScalarField rhoPCMmixPhi
(
    IOobject
    (
        "rhoPCMmix*phi",
        runTime.timeName(),
        mesh
    ),
    phi*fvc::interpolate(rhoPCMmix)
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
        (1.0 - beta*(T - T_S))
);

volScalarField DarcyTerm
(
    IOobject
    (
        "DarcyTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    DarcyConstant*pow(scalar(1)-gamma_L,2)/(pow(gamma_L,3)+DarcyTerm_epsilon)
);
//DarcyTerm.write();




#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"



volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
