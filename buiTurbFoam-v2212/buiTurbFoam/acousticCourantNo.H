scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

surfaceScalarField cSf(fvc::interpolate(Foam::sqrt(thermo.gamma() * thermo.p() / thermo.rho())) * mesh.magSf());
surfaceScalarField phiv((fvc::interpolate(rhoU) / fvc::interpolate(rho)) & mesh.Sf());
surfaceScalarField amaxSf(0.5 * (mag(phiv) + cSf));

if (mesh.nInternalFaces())
{
    scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum = 0.5 * gMax(sumAmaxSf / mesh.V().field()) * runTime.deltaTValue();

    meanCoNum = 0.5 * (gSum(sumAmaxSf) / gSum(mesh.V().field())) * runTime.deltaTValue();
}

Info << "Mean and max acoustic Courant Numbers = "
     << meanCoNum << " " << CoNum << endl;
