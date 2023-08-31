// Read field bounds
Info << "\nReading field bounds:" << endl;

// Pressure bounds
const dimensionedScalar pMin("pMin", p.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("pMin")));
Info << "pMin = " << pMin.value() << endl;

// Temperature bounds
const dimensionedScalar TMin("TMin", T.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("TMin")));
Info << "TMin = " << TMin.value() << endl;

// Enthalpy bounds
// const dimensionedScalar heMin("heMin", min(thermo.he((0 * thermo.p()) + pMin, (0 * thermo.T()) + TMin)));
// Info << "heMin = " << heMin.value() << endl;

// Density bounds
const dimensionedScalar rhoMin("rhoMin", rho.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("rhoMin")));
Info << "rhoMin = " << rhoMin.value() << endl;

Info << endl;
