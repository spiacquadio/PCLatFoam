IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reading density rho
dimensionedScalar rhoPCM(transportProperties.lookup("rhoPCM"));


// Reading liquid density
dimensionedScalar rhoeff_L(transportProperties.lookup("rhoeff_L"));
dimensionedScalar rhoPCM_L(transportProperties.lookup("rhoPCM_L"));

// Reading solid density
dimensionedScalar rhoeff_S(transportProperties.lookup("rhoeff_S"));
dimensionedScalar rhoPCM_S(transportProperties.lookup("rhoPCM_S"));

//Porosity
dimensionedScalar Porosity(transportProperties.lookup("Porosity"));

// Reading specific heat capacity
dimensionedScalar cpeff_S(transportProperties.lookup("cpeff_S"));
dimensionedScalar cpPCM_S(transportProperties.lookup("cpPCM_S"));

dimensionedScalar cpeff_L(transportProperties.lookup("cpeff_L"));
dimensionedScalar cpPCM_L(transportProperties.lookup("cpPCM_L"));

// Reading thermal conductivity
dimensionedTensor k_S(transportProperties.lookup("keff_S"));
dimensionedTensor k_L(transportProperties.lookup("keff_L"));

// Latent heat of melting
dimensionedScalar h_m(transportProperties.lookup("h_m"));

// Solidus and liquidus temperatures
dimensionedScalar T_S(transportProperties.lookup("T_S"));
dimensionedScalar T_L(transportProperties.lookup("T_L"));

// Calculate melting temperature
dimensionedScalar T_m
(
    "T_m",
    (T_L+T_S)/scalar(2.0)
);

// Reading parameters for Darcy source term
dimensionedScalar DarcyConstant(transportProperties.lookup("DarcyConstant"));
dimensionedScalar DarcyTerm_epsilon(transportProperties.lookup("DarcyTerm_epsilon"));

// Reading kinematic viscosity
dimensionedScalar nua(transportProperties.lookup("nua"));
dimensionedScalar nun(transportProperties.lookup("nun"));

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Reading Permeability [m²]
dimensionedTensor Permeability(transportProperties.lookup("Permeability"));
