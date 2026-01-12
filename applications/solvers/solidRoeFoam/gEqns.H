
rhsLm = fvc::surfaceIntegrate(flux_lm*magSf);

// if (angularMomentumConservation == "yes")
// {
//     rhsAm = fvc::surfaceIntegrate((xF ^ flux_lm)*magSf);
//     am.AMconservation(rhsLm, rhsLm1, rhsAm, stage);
// }

// Update coordinates
x += deltaT*(lm/rho);
xF += deltaT*(lmC/rho);
xN += deltaT*(lmN/rho);

lm += deltaT*rhsLm ;
lm -= deltaT*dampingCoeff*lm ;

// F += deltaT*fvc::surfaceIntegrate(flux_F*magSf);
F += deltaT*fvc::surfaceIntegrate((lmC/rho)*Sf);
