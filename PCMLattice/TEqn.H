{

    int iter = 0;

    scalar residual = 1;

    do
    {
        iter++;
        gamma_L.storePrevIter();

        // mixture density
        rhoeffmix=gamma_L*(rhoeff_L)
                +(scalar(1)-gamma_L)*(rhoeff_S);

        rhoPCMmix=gamma_L*(rhoPCM_L)
                +(scalar(1)-gamma_L)*(rhoPCM_S);

        // mixture specific heat
        cpeff = gamma_L*(cpeff_L) 
            + (scalar(1) - gamma_L)*(cpeff_S);

        cpPCM = gamma_L*(cpPCM_L) 
            + (scalar(1) - gamma_L)*(cpPCM_S);

        rhoPCMmixPhiCpPCM = fvc::interpolate(cpPCM)*phi*fvc::interpolate(rhoPCMmix);

        rhoeffmixPhi = fvc::interpolate(rhoeffmix)*phi;
        rhoPCMmixPhi = fvc::interpolate(rhoPCMmix)*phi;

        rhoeffmixCpeff = rhoeffmix*cpeff;

        // mixture thermal conductivity
        k = gamma_L*(k_L)
            + (scalar(1) - gamma_L)*(k_S);

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoeffmixCpeff,T)
          + fvm::div(rhoPCMmixPhiCpPCM, T)
          + h_m*Porosity*fvc::ddt(rhoPCMmix,gamma_L)
          + h_m*Porosity*fvc::div(rhoPCMmixPhi, gamma_L)
          - T_m*fvc::ddt(rhoeffmixCpeff)
          - T_m*fvc::div(rhoPCMmixPhiCpPCM)
          - fvm::laplacian(k, T)
        );

        TEqn.relax();

        TEqn.solve();

        Tcons = (T_L-T_S)*gamma_L+T_S;
        gamma_L = max(min(gamma_L + relaxGamma_L*cpeff/h_m*(T - Tcons),scalar(1)),scalar(0));

        residual = gMax(mag(gamma_L.primitiveField() - gamma_L.prevIter().primitiveField()));

        Info << "Correcting gamma_L, max residual = " << residual
        << endl;

    }while ((iter < minTiter || residual > gamma_L_Tol ) && iter <= maxTiter);

  
    rhok = (1.0 - beta*(T - T_S));
    nu = nua+(scalar(1)-gamma_L)*nun; //vel 0
    DarcyTerm = DarcyConstant*pow(scalar(1)-gamma_L,2)/(pow(gamma_L,3)+DarcyTerm_epsilon);


}
