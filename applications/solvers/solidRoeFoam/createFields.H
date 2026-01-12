Info << "Creating mesh parameters ..." << endl;

// Operations class
operations op(mesh);

// Point mesh
pointMesh pMesh(mesh);

// Material face area 
const surfaceScalarField& magSf = mesh.magSf();



// Material face area normal vector
const surfaceVectorField& Sf = mesh.Sf();
// Minimum edge length
const dimensionedScalar& h = op.minimumEdgeLength();

// Boundary patches
const polyBoundaryMesh& bm = mesh.boundaryMesh();
const label& symmetricPatchID = bm.findPatchID("symmetric");
const label& symmetricXpatchID = bm.findPatchID("symmetricX");
const label& symmetricYpatchID = bm.findPatchID("symmetricY");
const label& symmetricZpatchID = bm.findPatchID("symmetricZ");

Info << "Creating mesh coordinate fields ..." << endl;

// Material cell center coordinates
const volVectorField& C = mesh.C();

// Spatial cell center coordinates
volVectorField x
(
    IOobject("x", mesh),
    C
);

// Spatial nodal coordinates
pointVectorField xN
(
    IOobject("xN", mesh),
    pMesh,
    dimensionedVector("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
);
xN.primitiveFieldRef() = mesh.points();

// Material nodal coordinates
pointVectorField XN = xN;

// Spatial face center coordinates
surfaceVectorField xF = mesh.Cf();


Info << "Creating mesh normal fields ..." << endl;

// Material normals
const surfaceVectorField N = (Sf/mesh.magSf()).ref();

// Spatial normals
surfaceVectorField n = N;

// surfaceVectorField T1 = N;
// surfaceVectorField T2 = N;

surfaceVectorField T1
(
    IOobject("T1", mesh),
    mesh,
    dimensionedVector("T1", dimensionSet(0,0,0,0,0,0,0), vector::zero)
);
surfaceVectorField T2
(
    IOobject("T2", mesh),
    mesh,
    dimensionedVector("T2", dimensionSet(0,0,0,0,0,0,0), vector::zero)
);

// forAll(N, facei) 
// {
//     // Initialize an arbitrary vector
//     vector arbitraryVector;

//     // Choose a non-parallel arbitrary vector based on the normal vector components

//         arbitraryVector = vector(0, 0, 1); // z-axis
//     if (mag(arbitraryVector & N[facei]) != 0)    // If n is nearly aligned with z-axis
//     {
//         arbitraryVector = vector(0, 1, 0);   // Use the y-axis instead
//     }
//         // Info << "arbitraryVector" <<arbitraryVector<<endl;    
//     // First tangential vector (T1)
//     T1[facei] = N[facei] ^ arbitraryVector;
//         // Info << "T1" << T1[facei]<<endl;
//     T1[facei] /= mag(T1[facei]); // Normalize T1

//     // Second tangential vector (T2)
//     T2[facei] = N[facei] ^ T1[facei];
    
//     // Info << "T2" << T2[facei]<<endl;
//     T2[facei] /= mag(T2[facei]); // Normalize T2


//     // Now T1 and T2 are orthogonal to N and each other
// }
// Info << "N " << N <<endl;
// Info << "T1 " << T1<<endl;
// Info << "T2 " << T2<<endl;


// surfaceVectorField T1 = N;
// surfaceVectorField T2 = N;

// Loop through all faces in the mesh
forAll(N, facei) 
{
    // // Initialize an arbitrary vector
    vector arbitraryVector(1, 0, 0); // Default arbitrary vector along x-axis

    // Ensure the arbitrary vector is not parallel to the face normal N[facei]
    if (mag(arbitraryVector & N[facei]) > 0.999) // If nearly aligned
    {
        arbitraryVector = vector(0, 1, 0); // Switch to y-axis
        if (mag(arbitraryVector & N[facei]) > 0.999) // If still nearly aligned
        {
            arbitraryVector = vector(0, 0, 1); // Use z-axis as a last resort
        }
    }

    // Step 3: Project arbitraryVector onto the plane tangent to N[facei]
    T1[facei] = arbitraryVector - (arbitraryVector & N[facei]) * N[facei];

    // Find an orthogonal vector to N[facei]
    // T1[facei] = op.findOrthogonal(N[facei]);

    // Normalize T1 to make it a unit vector
    scalar T1Mag = mag(T1[facei]);
    if (T1Mag > SMALL) // Avoid division by zero
    {
        T1[facei] /= T1Mag;
    }
    else
    {
        FatalErrorInFunction << "Zero-length T1 encountered at face " << facei << abort(FatalError);
    }

    // Second tangential vector (T2): Cross product of N and T1
    T2[facei] = N[facei] ^ T1[facei];

    // Normalize T2 to make it a unit vector
    scalar T2Mag = mag(T2[facei]);
    if (T2Mag > SMALL) // Avoid division by zero
    {
        T2[facei] /= T2Mag;
    }
    else
    {
        FatalErrorInFunction << "Zero-length T2 encountered at face " << facei << abort(FatalError);
    }

    // Verify orthogonality of T1, T2, and N
    scalar tolerance = 1e-6;

    if (mag(N[facei] & T1[facei]) > tolerance) 
    {
        FatalErrorInFunction << "N and T1 are not orthogonal at face " << facei 
                             << ". Dot product: " << (N[facei] & T1[facei]) << abort(FatalError);
    }

    if (mag(N[facei] & T2[facei]) > tolerance) 
    {
        FatalErrorInFunction << "N and T2 are not orthogonal at face " << facei 
                             << ". Dot product: " << (N[facei] & T2[facei]) << abort(FatalError);
    }

    if (mag(T1[facei] & T2[facei]) > tolerance) 
    {
        FatalErrorInFunction << "T1 and T2 are not orthogonal at face " << facei 
                             << ". Dot product: " << (T1[facei] & T2[facei]) << abort(FatalError);
    }

    // Now T1 and T2 are guaranteed orthogonal to N and each other
}



// Info << "N " << N <<endl;
// Info << "T1 " << T1 <<endl;
// Info << "T2 " << T2 <<endl;


Info << "Creating linear momentum fields ..." << endl;

// Cell linear momentum
volVectorField lm
(
    IOobject
    (
        "lm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("lm", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
);

// Nodal linear momentum
pointVectorField lmN
(
    IOobject
    (
        "lmN",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh
);


Info << "Creating strain measure fields ..." << endl;

// Deformation gradient tensor
volTensorField F
(
    IOobject("F", mesh),
    mesh,
    tensor::I
);

//decompoz of tensor F
volVectorField Fx = op.decomposeTensorX(F);
volVectorField Fy = op.decomposeTensorY(F);
volVectorField Fz = op.decomposeTensorZ(F);

// Cofactor of deformation
volTensorField H
(
    IOobject("H", mesh),
    det(F)*op.invT(F)
);

// Jacobian of deformation
volScalarField J
(
    IOobject("J", mesh),
    det(F)
);



Info << "Creating constitutive model ..." << endl;

// Solid model class
solidModel model(F, mechanicalProperties);
model.correct();

// Density
const dimensionedScalar& rho = model.density();

// Pressure
volScalarField p = model.pressure();

// First Piola Kirchhoff stress tensor
volTensorField P = model.piola();
P.write();
volVectorField Px = op.decomposeTensorX(P);
volVectorField Py = op.decomposeTensorY(P);
volVectorField Pz = op.decomposeTensorZ(P);



volSymmTensorField sigma = symm((1.0/J)*(P & F.T()));

volScalarField sigmaEq
(
    "sigmaEq", sqrt((3.0/2.0)*magSqr(dev(sigma)))
);


Info << "Creating fields for wave speeds ..." << endl;

// Continuum mechanics class
mechanics mech(F, controlDict);

// Longitudinal wave speed
volScalarField Up
(
    IOobject("Up", mesh),
    mesh,
    model.Up()/beta           //?what is this beta? it's compressibility
);

// Shear wave speed
volScalarField Us
(
    IOobject ("Us", mesh),
    mesh,
    model.Us()*beta
);

mech.correct(Up, Us, F);

// Wave speed for time increment
volScalarField Up_time = (model.Up()/mech.stretch()).ref();

Info << "Creating fields for gradient ..." << endl;

// Gradient class
gradientSchemes grad(mesh,runTime);


// Gradient of cell linear momentum
volTensorField lmGrad = grad.gradient(lm);
// volTensorField lmGrad = fvc::grad(lm);

// Gradient of first Piola Kirchhoff stress tensor
volTensorField PxGrad = grad.gradient(Px);
volTensorField PyGrad = grad.gradient(Py);
volTensorField PzGrad = grad.gradient(Pz);
// compute F gradient
volTensorField FxGrad = grad.gradient(Fx);
volTensorField FyGrad = grad.gradient(Fy);
volTensorField FzGrad = grad.gradient(Fz);

// volTensorField PxGrad = fvc::grad(Px);
// volTensorField PyGrad = fvc::grad(Py);
// volTensorField PzGrad = fvc::grad(Pz);

// volTensorField FxGrad = fvc::grad(Fx);
// volTensorField FyGrad = fvc::grad(Fy);
// volTensorField FzGrad = fvc::grad(Fz);


Info << "Creating fields for reconstruction ..." << endl;
//limiter fields
volVectorField phi_lm
(
    IOobject("phi_lm", mesh),
    mesh,
    dimensionedVector("phi_lm", dimensionSet(0,0,0,0,0,0,0), vector::zero)
);

volTensorField phi_P
(
    IOobject("phi_P", mesh),
    mesh,
    dimensionedTensor("phi_P", dimensionSet(0,0,0,0,0,0,0), tensor::zero)
);
volTensorField phi_F
(
    IOobject("phi_F", mesh),
    mesh,
    dimensionedTensor("phi_F", dimensionSet(0,0,0,0,0,0,0), tensor::zero)
);

// Reconstruction of linear momentum
surfaceVectorField lm_M
(
    IOobject("lm_M", mesh),
    mesh,
    dimensionedVector("lm_M", lm.dimensions(), vector::zero)
);
surfaceVectorField lm_P = lm_M;
// grad.reconstruct(lm, lmGrad, lm_M, lm_P);


// Reconstruction of F tensor
surfaceTensorField F_M
(
    IOobject("F_M", mesh),
    mesh,
    tensor::I
);
surfaceTensorField F_P = F_M;
// grad.reconstruct(F, FxGrad, FyGrad, FzGrad, F_M, F_P);

// Reconstruction of PK1 stresses
surfaceTensorField P_M
(
    IOobject("P_M", mesh),
    mesh,
    dimensionedTensor("P_M", P.dimensions(), tensor::zero)
);
surfaceTensorField P_P = P_M;
// grad.reconstruct(P, PxGrad, PyGrad, PzGrad, P_M, P_P);

grad.reconstruct(lm, lmGrad, lm_M, lm_P, phi_lm);
grad.reconstruct(P, PxGrad, PyGrad, PzGrad, P_M, P_P, phi_P);
grad.reconstruct(F, FxGrad, FyGrad, FzGrad, F_M, F_P, phi_F);

// Reconstruction of traction
surfaceVectorField t_M
(
    IOobject("t_M", mesh),
    P_M & N
);
surfaceVectorField t_P
(
    IOobject("t_M", mesh),
    P_P & N
);
// surfaceVectorField t_P = (P_P & N).ref();

 //boundary firelds
Info << "Creating fields for riemann solver ..." << endl;

surfaceTensorField S_lm = mech.Smatrix_lm();
surfaceTensorField S_t = mech.Smatrix_t();
// Contact traction
surfaceVectorField tC = (0.5*(t_M + t_P) + (0.5*S_lm & (lm_P - lm_M))).ref();


volVectorField t_b
(
    IOobject
    (
        "t_b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField lm_b
(
    IOobject
    (
        "lm_b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Riemann solver
// #include "riemannSolver.H"



Info << "Creating fields for angular momentum ..." << endl;

// Angular momentum class
angularMomentum am(mesh, mechanicalProperties);

// RHS of linear momentum equation
volVectorField rhsLm
(
    IOobject("rhsLm", mesh),
    mesh,
    dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);

// RHS of linear momentum equation (stage 1)
volVectorField rhsLm1 = rhsLm;

// RHS of angular momentum equation
volVectorField rhsAm
(
    IOobject("rhsAm", mesh),
    mesh,
    dimensionedVector("rhsAm", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);

Info << "Creating fields for post-processing ..." << endl;

// Nodal displacements
pointVectorField uN
(
    IOobject
    (
        "uN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("uN", dimLength, vector::zero)
);
uN.write();

volVectorField u
(
    IOobject
    (
        "u",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("u", dimLength, vector::zero)
);
u.write();

volTensorField deformation
(
    IOobject
    (
        "deformation",
        runTime.timeName(), 
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    tensor::I
);
deformation.write();

Info << "Creating variables for time ..." << endl;

// Time increment
dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());

// Runge-Kutta stage
scalarList RKstages(2);
RKstages[0] = 0;
RKstages[1] = 1;


Info << "Printing data ..." << endl;

// Print material properties
model.printMaterialProperties();


// Print global linear and angular momentum
am.printGlobalMomentum(lm,x);

// Print centroid of geometry
mech.printCentroid();


//-------------Roe solver----------------

// Rotational matrix

surfaceTensorField R
(
    IOobject("R", mesh),
    mesh,
    tensor::I
);


surfaceTensorField RInvers
(
    IOobject("RInvers", mesh),
    mesh,
    tensor::I
);


surfaceTensorField RTranspos
(
    IOobject("RTranspos", mesh),
    mesh,
    tensor::I
);

surfaceTensorField RR
(
    IOobject("RR", mesh),
    mesh,
    tensor::I
);

surfaceTensorField RInversTranspos
(
    IOobject("RInversTranspos", mesh),
    mesh,
    tensor::I
);

surfaceTensorField RTransposInvers
(
    IOobject("RTransposInvers", mesh),
    mesh,
    tensor::I
);


//these are constant values
double epsilon = 1e-10; 

    // forAll(R, facei)
    // {
    //     R[facei].xx() = N[facei].x();
    //     R[facei].xy() = N[facei].y();
    //     R[facei].xz() = N[facei].z();
    //     R[facei].yx() = -(N[facei].y())/(std::sqrt(1-std::pow(N[facei].z(),2))+epsilon);
    //     R[facei].yy() =  (N[facei].x())/(std::sqrt(1-std::pow(N[facei].z(),2))+epsilon);
    //     R[facei].yz() = 0;
    //     R[facei].zx() = -(N[facei].x()*N[facei].z())/(std::sqrt(1-std::pow(N[facei].z(),2))+epsilon);
    //     R[facei].zy() = -(N[facei].y()*N[facei].z())/(std::sqrt(1-std::pow(N[facei].z(),2))+epsilon);
    //     R[facei].zz() = std::sqrt(1-std::pow(N[facei].z(),2));
    // }

        forAll(R, facei)
    {
        R[facei].xx() = N[facei].x();
        R[facei].xy() = N[facei].y();
        R[facei].xz() = N[facei].z();
        R[facei].yx() = T1[facei].x();
        R[facei].yy() = T1[facei].y();
        R[facei].yz() = T1[facei].z();
        R[facei].zx() = T2[facei].x();
        R[facei].zy() = T2[facei].y();
        R[facei].zz() = T2[facei].z();
    }


    //     forAll(R, facei)
    // {

    //     R[facei].xx() = N[facei].x();
    //     R[facei].xy() = T1[facei].x();
    //     R[facei].xz() = T2[facei].x();
    //     R[facei].yx() = N[facei].y();
    //     R[facei].yy() = T1[facei].y();
    //     R[facei].yz() = T2[facei].y(); 
    //     R[facei].zx() = N[facei].z();
    //     R[facei].zy() = T1[facei].z();
    //     R[facei].zz() = T2[facei].z();
    // }





// RInvers = inv(R);
RTranspos = R.T();
// RInversTranspos = RInvers.T();
// RTransposInvers = inv(RTranspos);



// Info << "R  " << R <<endl;
// const surfaceVectorField N_Rotated = R & N ;
// const surfaceVectorField N_Rotated_back = RTranspos  & N_Rotated;
// Info << "N  " << N <<endl;
// Info << "N_Rotated  " << N_Rotated <<endl;
// Info << "N_Rotated_back  " << N_Rotated_back <<endl;
// Info << "N-N_Rotated_back  " << N-N_Rotated_back <<endl;
// // Info << "RTranspos  " << RTranspos <<endl;
// // Info << "RInvers  " << RInvers <<endl;
// // Info << "RInversTranspos  " << RInversTranspos <<endl;
// Info << "RR  " << RR <<endl;
// Info << "RR_Rotated  " << RR_Rotated <<endl;
// Info << "RR_Rotated_back " << RR_Rotated_back <<endl;
// Info << "RR-RR_Rotated_back = " << RR-RR_Rotated_back <<endl;

//--------------------------

//apply rotationl inveriant to the reconstructed fields
// surfaceVectorField lm_M_hat = (R & lm_M).ref();
// surfaceVectorField lm_P_hat = (R & lm_P).ref();

// surfaceTensorField F_M_hat = (R & F_M & RTranspos).ref();
// surfaceTensorField F_P_hat = (R & F_P & RTranspos).ref();
surfaceVectorField lm_M_hat = (R & lm_M);
surfaceVectorField lm_P_hat = (R & lm_P);

surfaceTensorField F_M_hat = (R & F_M & RTranspos);
surfaceTensorField F_P_hat = (R & F_P & RTranspos);

// Reconstruction of PK1 stresses
//! first we have to define PK1 stress in the roated fram
// Solid model class
// solidModelPlusMinus modelPlusMinus (F_P_hat, F_M_hat ,mechanicalProperties);
// modelPlusMinus.correct();

// surfaceTensorField P_hat_P = (modelPlusMinus.piola_P());
// surfaceTensorField P_hat_M = (modelPlusMinus.piola_M());

//? what is instead I rotated the Pk1 instead ? lik
// surfaceTensorField P_hat_P = (R & P_P & RTranspos).ref();
// surfaceTensorField P_hat_M = (R & P_M & RTranspos).ref();

surfaceTensorField P_hat_P = (R & P_P & RTranspos);
surfaceTensorField P_hat_M = (R & P_M & RTranspos);


// //---------------------------
// define the fluxes after roation on the left and right  
surfaceVectorField flux_lm_M_hat
(
    IOobject("flux_lm_M_hat", mesh),
    mesh,
    dimensionedVector("flux_lm_M_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);

surfaceVectorField flux_lm_P_hat
(
    IOobject("flux_lm_P_hat", mesh),
    mesh,
    dimensionedVector("flux_lm_P_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)

);

surfaceTensorField flux_F_M_hat
(
    IOobject("flux_F_M_hat", mesh),
    mesh,
    dimensionedTensor("flux_F_M_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)

);
surfaceTensorField flux_F_P_hat
(
    IOobject("flux_F_P_hat", mesh),
    mesh,
    dimensionedTensor("flux_F_P_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)

);

//------assign the flux value according to the rotational inveriant----
//-------------------------------------------------------------------------------------
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
    


//-------Roe avarage values-----------
//difference fields

// surfaceVectorField deltaF_lm = (lm_P_hat - lm_M_hat).ref();
// surfaceTensorField deltaF_F  = (F_P_hat - F_M_hat).ref();
surfaceVectorField deltaF_lm = (lm_P_hat - lm_M_hat);
surfaceTensorField deltaF_F  = (F_P_hat - F_M_hat);
  
//rho_ is defined at the surface for ease access
surfaceScalarField rho_
(
    IOobject("rho_", mesh),
    mesh,
    model.density()
);

//Eigenvalues
surfaceScalarField Up_
(
    IOobject("Up_", mesh),
    mesh,
    model.Up()
);
// Info << "Dimensions of Up_: " << Up_.dimensions() << endl;


surfaceScalarField Us_
(
    IOobject("Us_", mesh),
    mesh,
    model.Us()
);

surfaceScalarField lambda_ //lame's 1st parameter
(
    IOobject("lambda_", mesh),
    mesh,
    model.lambda()
);
// Info << "Dimensions of Us_: " << Us_.dimensions() << endl;

//eigenvector linear momentum
surfaceVectorField l1_lm
(
    IOobject("l1_lm", mesh),
    mesh,
    vector::zero
);

surfaceVectorField l2_lm
(
    IOobject("l2_lm", mesh),
    mesh,
    vector::zero
);

surfaceVectorField l3_lm
(
    IOobject("l3_lm", mesh),
    mesh,
    vector::zero
);
surfaceVectorField l4_lm
(
    IOobject("l4_lm", mesh),
    mesh,
    vector::zero
);

surfaceVectorField l5_lm
(
    IOobject("l5_lm", mesh),
    mesh,
    vector::zero
);

surfaceVectorField l6_lm
(
    IOobject("l6_lm", mesh),
    mesh,
    vector::zero
);

//eigenvector deformation tensor presented in tensor form
surfaceTensorField l1_F
(
    IOobject("l1_F", mesh),
    mesh,
    dimensionedTensor("l1_F", dimensionSet(-1,2,1,0,0,0,0), tensor::zero)
);

surfaceTensorField l2_F
(
    IOobject("l2_F", mesh),
    mesh,
    dimensionedTensor("l1_F", dimensionSet(-1,2,1,0,0,0,0), tensor::zero)
);

surfaceTensorField l3_F
(
    IOobject("l3_F", mesh),
    mesh,
    dimensionedTensor("l1_F", dimensionSet(-1,2,1,0,0,0,0), tensor::zero)
);
surfaceTensorField l4_F
(
    IOobject("l4_F", mesh),
    mesh,
    dimensionedTensor("l1_F", dimensionSet(-1,2,1,0,0,0,0), tensor::zero)
);

surfaceTensorField l5_F
(
    IOobject("l5_F", mesh),
    mesh,
    dimensionedTensor("l1_F", dimensionSet(-1,2,1,0,0,0,0), tensor::zero)
);

surfaceTensorField l6_F
(
    IOobject("l6_F", mesh),
    mesh,
    dimensionedTensor("l1_F", dimensionSet(-1,2,1,0,0,0,0), tensor::zero)
);

//!these are constants all the time
forAll(l1_lm, faceI)
{
    l1_lm[faceI][0] = 1.0; // x-component

    l2_lm[faceI][0] = 1.0; // x-component

    l3_lm[faceI][1] = 1.0; // Y-component

    l4_lm[faceI][2] = 1.0; // Z-component

    l5_lm[faceI][1] = 1.0; // Y-component

    l6_lm[faceI][2] = 1.0; // Z-component


    //assign values to 
     l1_F[faceI][0] = -1/( rho_[faceI] * Up_[faceI]); // xx-component
     l2_F[faceI][0] =  1/( rho_[faceI] * Up_[faceI]); // xx-component

     l3_F[faceI][3] = -1/( rho_[faceI] * Us_[faceI]); // zx-component
     l4_F[faceI][6] = -1/( rho_[faceI] * Us_[faceI]); // yx-component

     l5_F[faceI][3] =  1/( rho_[faceI] * Us_[faceI]); // zx-component
     l6_F[faceI][6] =  1/( rho_[faceI] * Us_[faceI]); // yx-component


}

//wave strength
surfaceScalarField alpha1
(
    IOobject
    ("alpha1",mesh),
    mesh,
    dimensionedScalar("alpha1", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);

surfaceScalarField alpha2
(
    IOobject
    ("alpha2",mesh),
    mesh,
    dimensionedScalar("alpha2", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);

surfaceScalarField alpha3
(
    IOobject
    ("alpha3",mesh),
    mesh,
    dimensionedScalar("alpha3", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);

surfaceScalarField alpha4
(
    IOobject
    ("alpha4",mesh),
    mesh,
    dimensionedScalar("alpha4", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);

surfaceScalarField alpha5
(
    IOobject
    ("alpha5",mesh),
    mesh,
    dimensionedScalar("alpha5", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);

surfaceScalarField alpha6
(
    IOobject
    ("alpha6",mesh),
    mesh,
    dimensionedScalar("alpha6", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);

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
surfaceVectorField sum_lm
(
    IOobject("sum_lm", mesh),
    mesh,
    dimensionedVector("flux_lm_M_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);
// sum_lm = ( l1_lm );
sum_lm = (Up_ * alpha1 * l1_lm )
       + (Up_ * alpha2 * l2_lm )
       + (Us_ * alpha3 * l3_lm )
       + (Us_ * alpha4 * l4_lm )
       + (Us_ * alpha5 * l5_lm )
       + (Us_ * alpha6 * l6_lm ) ;

// Info << sum_lm <<endl;

surfaceTensorField sum_F
(
    IOobject("sum_F", mesh),
    mesh,
    dimensionedTensor("flux_F_P_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
);

sum_F  = (Up_ * alpha1 * l1_F )
       + (Up_ * alpha2 * l2_F )
       + (Us_ * alpha3 * l3_F )
       + (Us_ * alpha4 * l4_F )
       + (Us_ * alpha5 * l5_F )
       + (Us_ * alpha6 * l6_F )  ;

//Roe Fluxs
surfaceVectorField roeFlux_lm_hat
(
    IOobject("roeFlux_lm_hat", mesh),
    mesh,
    dimensionedVector("flux_lm_M_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);

surfaceTensorField roeFlux_F_hat
(
    IOobject("roeFlux_F_hat", mesh),
    mesh,
    dimensionedTensor("flux_F_M_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
);

// roeFlux_lm_hat = 0.5*  (flux_lm_P_hat + flux_lm_M_hat  - sum_lm);
// roeFlux_F_hat  = 0.5*  (flux_F_P_hat  + flux_F_M_hat   - sum_F);


roeFlux_lm_hat = 0.5*  (flux_lm_P_hat + flux_lm_M_hat  + sum_lm);
roeFlux_F_hat  = 0.5*  (flux_F_P_hat  + flux_F_M_hat   + sum_F);

// Info<< "roeFlux_lm_hat" <<roeFlux_lm_hat<<endl;

//rotate the system back
surfaceVectorField  flux_lm = RTranspos & roeFlux_lm_hat;
// surfaceTensorField flux_F = RInvers & roeFlux_F_hat & RInversTranspos;
surfaceTensorField  flux_F = (RTranspos & roeFlux_F_hat) & R;


//apply boundary conditions
// in contact algorithm. the bcs are applied to tc and lmc
// Compute boundary values
lm_b.correctBoundaryConditions();
t_b.correctBoundaryConditions();

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


// #include "riemannSolver.H"

//apply the curl free correction algorithm
surfaceScalarField N_norm_squared = sqr(mag(N));

surfaceVectorField lmC
(
    IOobject("lmC", mesh),
    mesh,
    dimensionedVector("lmC", dimensionSet(1 ,-2 ,-1 ,0 ,0 ,0 ,0), vector::zero)
);

lmC = (rho_/N_norm_squared) * (flux_F & N) ;

Info << "Creating fields for the constrained procedure ..." << endl;

// Constrained class
interpolationSchemes interpolate(mesh);

// Cell averaged linear momentum
volVectorField lmR = interpolate.surfaceToVol(lmC);

// Local gradient of cell averaged linear momentum
volTensorField lmRgrad = grad.localGradient(lmR, lmC);

// Constrained fluxes
interpolate.volToPoint(lmR, lmRgrad, lmN);
#include "strongBCs.H"
lmN.correctBoundaryConditions();
lmC = interpolate.pointToSurface(lmN);

//updating the flux after correction:
flux_F = (1/rho_)* (lmC * N);

