volScalarField &p = thermo.p();
volScalarField &T = thermo.T();
const volScalarField &psi = thermo.psi();
const volScalarField &mu = thermo.mu();

bool inviscid(true);
if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}
