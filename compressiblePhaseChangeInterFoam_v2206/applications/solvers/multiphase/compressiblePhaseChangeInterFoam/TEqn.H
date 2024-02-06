{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             divUp()// - contErr/rho*p
          +  (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))() - contErr*K
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
     ==
        fvOptions(rho, T)
        - mEvaporation*vaporizationEnthalpy/(alpha1*mixture.thermo1().Cv()()+alpha2*mixture.thermo2().Cv()()) //*mixture.thermo1().rho()*cellVolu
        + mCondensation*vaporizationEnthalpy/(alpha1*mixture.thermo1().Cv()()+alpha2*mixture.thermo2().Cv()()) // *mixture.thermo2().rho()*cellVolu
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correctThermo();
    mixture.correct();
}
