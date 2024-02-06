#include "createRDeltaT.H"

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

Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);


dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt(alpha1*fvc::div(phi));

#include "createAlphaFluxes.H"

// Construct compressible turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi10,
    mixture
);

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"



Info<< "Reading thermalConductivity\n" << endl;

IOdictionary thermalConductivity
(
    IOobject
    (
        "thermalConductivity",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading field kLaminar\n" << endl;
volScalarField kLaminar
(
    IOobject
    (
        "kLaminar",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroKLaminar", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0)
);

Info<< "Reading field Tsaturation\n" << endl;
volScalarField Tsaturation
(
    IOobject
    (
        "Tsaturation",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroT", dimTemperature, 0.0)
);

// adimensionalize unit for pressure
dimensionedScalar PaUnit
(
    "PaUnit",
    dimPressure,
    1
);

// adimensionalize unit for Temperature
dimensionedScalar KelvinUnit
(
    "KelvinUnit",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    1
);


// Phase 1 thermal conductivity
dimensionedScalar k1
(
    "k1",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    thermalConductivity.lookup("k1")
);

// Phase 2 thermal conductivity
dimensionedScalar k2
(
    "k2",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    thermalConductivity.lookup("k2")
);


Info<< "Reading phaseChangeProperties\n" << endl;

IOdictionary phaseChangeProp
(
    IOobject
    (
        "phaseChangeProp",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar betaEvaporation
(
    "betaEvaporation",
    dimensionSet(0, 0, -1, 0, 0, 0, 0),
    phaseChangeProp.lookup("betaEvaporation")
);

dimensionedScalar betaCondensation
(
    "betaCondensation",
    dimensionSet(0, 0, -1, 0, 0, 0, 0),
    phaseChangeProp.lookup("betaCondensation")
);

dimensionedScalar Tsat
(
    "Tsat",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    phaseChangeProp.lookup("Tsat")
);

dimensionedScalar vaporizationEnthalpy
(
    "vaporizationEnthalpy",
    dimensionSet(0, 2, -2, 0, 0, 0, 0),
    phaseChangeProp.lookup("vaporizationEnthalpy")
);

Info<< "Reading field mEvaporation\n" << endl;
volScalarField mEvaporation
(
    IOobject
    (
        "mEvaporation",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroMCondensation", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
);

Info<< "Reading field mCondensation\n" << endl;
volScalarField mCondensation
(
    IOobject
    (
        "mCondensation",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroMCondensation", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
);



volScalarField cellVolu
(
IOobject
    (
    "cellVolu",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("cellVolu", dimVolume, 0.0)
);   
cellVolu.ref() = mesh.V(); 

Info<< "Reading field kGradT\n" << endl;
volVectorField kGradT
(
    IOobject
    (
        "kGradT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("kGradTzero", dimensionSet(1, 0, -3, 0, 0, 0, 0), Foam::vector(0,0,0))
);
