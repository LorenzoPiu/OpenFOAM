tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

Tsaturation = (273.25 - 124.8 + 91.08*pow(p/100000/PaUnit,(1/4.629)))*KelvinUnit;
mEvaporation = betaEvaporation*alpha1*mixture.thermo1().rho()*(max(T-Tsaturation, T*0)/Tsaturation);
mCondensation = betaCondensation*alpha2*mixture.thermo1().rho()*(max(Tsaturation-T, T*0)/Tsaturation);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();

    talphaPhi1 = new surfaceScalarField
    (
        IOobject
        (
            "alphaPhi1",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(alphaPhi10.dimensions(), Zero)
    );

    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;

const surfaceScalarField& alphaPhi1 = talphaPhi1();
surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

volScalarField::Internal contErr
(
    (
        fvc::ddt(rho) + fvc::div(rhoPhi)
      - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
      - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
    )()
);
