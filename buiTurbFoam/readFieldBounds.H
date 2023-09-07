// Read field bounds
Info << "\nReading field bounds:" << endl;

// Pressure bounds
const dimensionedScalar pMin("pMin", p.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("pMin")));
const dimensionedScalar pMax("pMax", p.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("pMax")));
Info << "pMin = " << pMin.value() << endl;
Info << "pMax = " << pMax.value() << endl;

// Temperature bounds
const dimensionedScalar TMin("TMin", T.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("TMin")));
const dimensionedScalar TMax("TMax", T.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("TMax")));
Info << "TMin = " << TMin.value() << endl;
Info << "TMax = " << TMax.value() << endl;

// Enthalpy bounds
const dimensionedScalar heMin("heMin", min(thermo.he((0 * thermo.p()) + pMin, (0 * thermo.T()) + TMin)));
const dimensionedScalar heMax("heMax", max(thermo.he((0 * thermo.p()) + pMax, (0 * thermo.T()) + TMax)));
Info << "heMin = " << heMin.value() << endl;
Info << "heMax = " << heMax.value() << endl;

// Density bounds
const dimensionedScalar rhoMin("rhoMin", rho.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("rhoMin")));
const dimensionedScalar rhoMax("rhoMax", rho.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("rhoMax")));
Info << "rhoMin = " << rhoMin.value() << endl;
Info << "rhoMax = " << rhoMax.value() << endl;

Info << endl;
