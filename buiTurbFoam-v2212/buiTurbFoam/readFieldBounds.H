// Read field bounds
Info << "Reading field bounds\n"
     << endl;

// Pressure bounds
const dimensionedScalar pMin("pMin", p.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("pMin")));
Info << "\n pMin = " << pMin.value() << endl;

// Temperature bounds
const dimensionedScalar TMin("TMin", T.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("TMin")));
Info << "\n TMin = " << TMin.value() << endl;

// Enthalpy bounds
const dimensionedScalar heMin("heMin", min(thermo.Cp()) * TMin);
Info << "\n heMin = " << heMin.value() << endl;

// Density bounds
const dimensionedScalar rhoMin("rhoMin", rho.dimensions(), readScalar(mesh.solutionDict().subDict("fieldBounds").lookup("rhoMin")));
Info << "\n rhoMin = " << rhoMin.value() << endl;
