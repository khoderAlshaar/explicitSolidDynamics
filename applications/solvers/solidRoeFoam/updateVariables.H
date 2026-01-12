// Strain measures
H = det(F)*op.invT(F);
J = det(F);

// Wave speeds
mech.correct(Up, Us, F);
Up_time = model.Up()/mech.stretch();

// Constitutive model
model.correct();
P = model.piola();

// Cell gradients
lmGrad = grad.gradient(lm);
grad.gradient(P, PxGrad, PyGrad, PzGrad);
grad.gradient(F, FxGrad, FyGrad, FzGrad);


// Reconstruction
// grad.reconstruct(lm, lmGrad, lm_M, lm_P);
// grad.reconstruct(P, PxGrad, PyGrad, PzGrad, P_M, P_P);
// grad.reconstruct(F, FxGrad, FyGrad, FzGrad, F_M, F_P);


grad.reconstruct(lm, lmGrad, lm_M, lm_P, phi_lm);
grad.reconstruct(P, PxGrad, PyGrad, PzGrad, P_M, P_P, phi_P);
grad.reconstruct(F, FxGrad, FyGrad, FzGrad, F_M, F_P, phi_F);

// this is being updated for the boundary conditions requiermant
t_M = P_M & N;
t_P = P_P & N;
S_lm = mech.Smatrix_lm();
S_t = mech.Smatrix_t();
// Riemann solver
#include "riemannSolverRoe.H" 

// Info << "lm flux (tC)" << flux_lm <<endl;

// Nodal linear momentum
 N_norm_squared = sqr(mag(N));

 lmC = (rho_/N_norm_squared) * (flux_F & N) ;
//  lmC = (rho_) * (flux_F & N) ;

//C-TOUCH algorithm
//step a
lmR = interpolate.surfaceToVol(lmC);
lmRgrad = grad.localGradient(lmR, lmC);
interpolate.volToPoint(lmR, lmRgrad, lmN);
#include "strongBCs.H"
lmN.correctBoundaryConditions();
//step d
// Constrained fluxes
lmC = interpolate.pointToSurface(lmN);

// Info << "lmC after correction" << lmC <<endl;

//updating the flux after correction:
flux_F = (1/rho_)* (lmC * N);


// Info << "flux_F = (1/rho_)* (lmC * N); updateVariables  " << flux_F - (1/rho_)* (lmC * N)<<endl;