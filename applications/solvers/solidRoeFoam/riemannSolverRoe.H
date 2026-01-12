// Roe Riemann solver
//apply rotationl inveriant to the reconstructed fields

 lm_M_hat = R & lm_M;
 lm_P_hat = R & lm_P;

//  Info<< "lm_M " <<lm_M<<endl;
//  Info<< "lm_P " <<lm_P<<endl;
//   Info<< "R " <<R<<endl;
//  Info<< "lm_P_hat " <<lm_P_hat<<endl;


 F_M_hat = (R & F_M) & RTranspos;
 F_P_hat = (R & F_P) & RTranspos;

//  Info <<"F_M "<<F_M<<endl;
//  Info <<"-------------------"<<endl;
//  Info <<"F_M_hat"<<F_M_hat<<endl;

// modelPlusMinus.correct();
// P_hat_P = modelPlusMinus.piola_P();
// P_hat_M = modelPlusMinus.piola_M();
 P_hat_M = (R & P_M) & RTranspos;
 P_hat_P = (R & P_P) & RTranspos;

 //difference fields
deltaF_lm  = lm_P_hat - lm_M_hat;
deltaF_F   = F_P_hat -  F_M_hat;

//linear momentum flux
    forAll(flux_lm_P_hat, facei)
    {
        flux_lm_P_hat[facei].x() = P_hat_P[facei].xx();
        flux_lm_P_hat[facei].y() = P_hat_P[facei].yx();
        flux_lm_P_hat[facei].z() = P_hat_P[facei].zx();
        
        flux_lm_M_hat[facei].x() = P_hat_M[facei].xx();
        flux_lm_M_hat[facei].y() = P_hat_M[facei].yx();
        flux_lm_M_hat[facei].z() = P_hat_M[facei].zx();

        //deformation gradient tensor flux
        flux_F_M_hat[facei].xx() = (1/rho.value()) * lm_M_hat[facei].x();
        flux_F_M_hat[facei].yx() = (1/rho.value()) * lm_M_hat[facei].y();
        flux_F_M_hat[facei].zx() = (1/rho.value()) * lm_M_hat[facei].z();
       
        flux_F_P_hat[facei].xx() = (1/rho.value()) * lm_P_hat[facei].x();
        flux_F_P_hat[facei].yx() = (1/rho.value()) * lm_P_hat[facei].y();
        flux_F_P_hat[facei].zx() = (1/rho.value()) * lm_P_hat[facei].z();
    }
// Info << flux_lm_M_hat <<endl;


 //eigenvalues and eigenvectors are constants no needd to update them again for a linear system


forAll(alpha1,facei)
{
    alpha1[facei] = 0.5*( - Up_[facei] * rho_[facei] * deltaF_F[facei].xx() - (lambda_[facei]/Up_[facei])*(deltaF_F[facei].yy()+ deltaF_F[facei].zz() )  +  deltaF_lm[facei].x() );
    alpha2[facei] = 0.5*(   Up_[facei] * rho_[facei] * deltaF_F[facei].xx() + (lambda_[facei]/Up_[facei])*(deltaF_F[facei].yy()+ deltaF_F[facei].zz() )  +  deltaF_lm[facei].x() );

    alpha3[facei] = 0.5*( - Us_[facei]* rho_[facei] * (deltaF_F[facei].xy() +deltaF_F[facei].yx() ) + deltaF_lm[facei].y() );
    alpha4[facei] = 0.5*( - Us_[facei]* rho_[facei] * (deltaF_F[facei].xz() +deltaF_F[facei].zx() ) + deltaF_lm[facei].z() );
    
    alpha5[facei] = 0.5*(   Us_[facei]* rho_[facei] * (deltaF_F[facei].xy() +deltaF_F[facei].yx() ) + deltaF_lm[facei].y() );
    alpha6[facei] = 0.5*(   Us_[facei]* rho_[facei] * (deltaF_F[facei].xz() +deltaF_F[facei].zx() ) + deltaF_lm[facei].z() );
    

}

//sum
sum_lm = (Up_ * alpha1 * l1_lm )
       + (Up_ * alpha2 * l2_lm )
       + (Us_ * alpha3 * l3_lm )
       + (Us_ * alpha4 * l4_lm )
       + (Us_ * alpha5 * l5_lm )
       + (Us_ * alpha6 * l6_lm ) ;

sum_F  = (Up_ * alpha1 * l1_F )
       + (Up_ * alpha2 * l2_F )
       + (Us_ * alpha3 * l3_F )
       + (Us_ * alpha4 * l4_F )
       + (Us_ * alpha5 * l5_F )
       + (Us_ * alpha6 * l6_F )  ;


roeFlux_lm_hat = 0.5*  (flux_lm_P_hat + flux_lm_M_hat  + sum_lm);
roeFlux_F_hat  = 0.5*  (flux_F_P_hat  + flux_F_M_hat   + sum_F);

// Info<< "flux_lm_P_hat " <<flux_lm_P_hat<<endl;
// Info<< "flux_lm_M_hat " <<flux_lm_M_hat<<endl;
// Info<< "sum_lm " <<flux_lm_M_hat<<endl;
// Info<< "roeFlux_lm_hat " <<roeFlux_lm_hat<<endl;

//rotate the system back
 flux_lm = RTranspos & roeFlux_lm_hat;
//  flux_F = (RInvers & roeFlux_F_hat) & RInversTranspos;
 flux_F = (RTranspos & roeFlux_F_hat) & R;


// Info <<" RInvers  "<<  RInvers<< endl;
// Info <<" RTranspos  "<<  RTranspos<< endl;
// Info <<" flux_lm "<<  flux_lm<< endl;
// Info <<" flux_F "<<  flux_F<< endl;

// Compute boundary values
lm_b.correctBoundaryConditions();
t_b.correctBoundaryConditions();

// Info << "lm_b  " <<lm_b<<endl;

forAll(mesh.boundary(), patchi)
{

        forAll(mesh.boundary()[patchi], facei)
        {
            flux_lm.boundaryFieldRef()[patchi][facei] =  
                t_b.boundaryField()[patchi][facei];

            flux_F.boundaryFieldRef()[patchi][facei] =  (1/rho_.boundaryFieldRef()[patchi][facei])*
                ( lm_b.boundaryField()[patchi][facei] * N.boundaryField()[patchi][facei]);
    
        }

}
// Info << "flux_F_Riemann solver  " <<flux_F<<endl;