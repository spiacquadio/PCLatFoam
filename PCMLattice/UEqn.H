    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
        + (1/Porosity)*fvm::div(phi, U)
        - fvm::laplacian(nu_eff, U)
        + fvm::Sp(DarcyTerm*Porosity/rhoPCM, U)
        + ((inv(Permeability))&(fvc::Sp(nu*Porosity, U)))
        + MRF.DDt(U)
        ==
        fvModels.source(U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - Porosity * ghf*fvc::snGrad(rhok)
                  - Porosity * fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
        
        fvConstraints.constrain(U);
    }
