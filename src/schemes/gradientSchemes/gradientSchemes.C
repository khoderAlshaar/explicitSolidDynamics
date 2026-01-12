/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gradientSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //


gradientSchemes::gradientSchemes
(
    const fvMesh& vm, const Time& runTime
)
:
    mesh_(vm),
    own_(mesh_.owner()),
    nei_(mesh_.neighbour()),
    X_(mesh_.C()),
    XF_(mesh_.Cf()),
    XN_(mesh_.points()),
 
    Ainv_
    (
        IOobject("Ainv", mesh_),
        mesh_,
        dimensionedTensor("Ainv", dimensionSet(0,2,0,0,0,0,0), tensor::zero)
    ),

    AinvLocal_
    (
        IOobject("AinvLocal", mesh_),
        mesh_,
        dimensionedTensor
        (
            "AinvLocal",
            dimensionSet(0,2,0,0,0,0,0),
            tensor::zero
        )
    )
    ,
    
    //     // Finite volume solution dictionary
    //reading dicts
    runParameters_
    (
        IOobject
        (
            "runParameters",
            runTime.constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
        order_
    (
        runParameters_.lookupOrAddDefault<word>("reconstructionOrder", "first")
        // runParameters_.lookup("reconstructionOrder")
    ),
        limiter_
    (
        runParameters_.lookupOrAddDefault<word>("limiter", "noLimiter")
        // runParameters_.lookup("limiter")
    )



{
    gradientSchemes::distanceMatrix(Ainv_);
    gradientSchemes::distanceMatrixLocal(AinvLocal_);


    //  order_ = fvSolution_.lookup("reconstructionOrder");
     if
    (
        order_ != "first" && order_ != "second"
    )
    {
        FatalErrorIn("gradientSchemes.H")
            << "Valid type entries are 'first' or 'second' "
            << "for reconstructionOrder"
            << abort(FatalError);
    }
     if
    (
        limiter_ != "no" && limiter_ != "yes"
    )
    {
        FatalErrorIn("gradientSchemes.H")
            << "Valid type entries are 'no' or 'yes' "
            << "for limiter"
            << abort(FatalError);
    }
    Info << " reconstruction scheme order: " << order_ <<endl;
    Info << " Limiter: " << limiter_ <<endl;

}

// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

gradientSchemes::~gradientSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradientSchemes::distanceMatrix
(
    GeometricField<tensor, fvPatchField, volMesh>& U
)
{
    forAll(own_, faceID)
    {
        const label& ownCellID = own_[faceID];
        const label& neiCellID = nei_[faceID];
        const vector& dOwn = X_[neiCellID] - X_[ownCellID];
        const vector& dNei  = X_[ownCellID] - X_[neiCellID];

        U[ownCellID] += dOwn*dOwn;
        U[neiCellID] += dNei*dNei;
    }

    if (Pstream::parRun())
    {
        forAll(mesh_.boundary(), patchID)
        {
            if (mesh_.boundary()[patchID].coupled())
            {
                tmp<vectorField> tmp_X_nei = X_.boundaryField()[patchID].patchNeighbourField();
                const vectorField& X_nei = tmp_X_nei();

                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    const vector& d = X_nei[facei] - X_[bCellID];
                    U[bCellID] += d*d;
                }
            }
        }

        U.correctBoundaryConditions();
    }

    U.primitiveFieldRef() = inv(U.internalField());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::distanceMatrixLocal
(
    GeometricField<tensor, fvPatchField, volMesh>& Ainv
) const
{
    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tvf
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject("distanceMatrixLocal", mesh_),
            mesh_,
            dimensioned<tensor>("0", Ainv.dimensions(), pTraits<tensor>::zero)
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> dCd = tvf();

    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];
        const vector& dOwn = XF_[faceID] - X_[ownID];
        const vector& dNei = XF_[faceID] - X_[neiID];

        dCd[ownID] += dOwn*dOwn;
        dCd[neiID] += dNei*dNei;
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundary()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            vector d = XF_.boundaryField()[patchID][facei] - X_[bCellID];
            dCd[bCellID] += d*d;

            //! works only with 3D geometry. Further invistigation required
            // if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            // {
            //     const label& faceID =
            //         mesh_.boundary()[patchID].start() + facei;

            //     forAll(mesh_.faces()[faceID], nodei)
            //     {
            //         const label& nodeID = mesh_.faces()[faceID][nodei];

            //         d = XN_[nodeID] - X_[bCellID];
            //         dCd[bCellID] += d * d;

            //         for (int i=0; i<7; i++)
            //         {
            //             d =
            //                 ((((i+1)*XN_[nodeID])
            //               + ((7 - i)*XF_.boundaryField()[patchID][facei]))/8.0)
            //               - X_[bCellID];
            //             dCd[bCellID] += d * d;
            //         }
            //     }
            // }
        }
    }

    Ainv.primitiveFieldRef() = inv(dCd.internalField());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volVectorField gradientSchemes::gradient
(
    const GeometricField<scalar, fvPatchField, volMesh>& U
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ugrad = tvf();

    forAll(own_, faceID)
    {
        const label& cellID = own_[faceID];
        const label& neiID = nei_[faceID];
        const vector& dcell = X_[neiID] - X_[cellID];
        const vector& dnei = X_[cellID] - X_[neiID];

        Ugrad[cellID] += Ainv_[cellID] & (U[neiID] - U[cellID])*dcell;
        Ugrad[neiID] += Ainv_[neiID] & (U[cellID] - U[neiID])*dnei;
    }

    if (Pstream::parRun())
    {
        forAll(mesh_.boundary(), patchID)
        {
            if (mesh_.boundary()[patchID].coupled())
            {
                tmp<vectorField> tmp_X_nei = X_.boundaryField()[patchID].patchNeighbourField();
                const vectorField& X_nei = tmp_X_nei();

                tmp<scalarField> tmp_U_nei = U.boundaryField()[patchID].patchNeighbourField();
                const scalarField& U_nei = tmp_U_nei();

                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    const vector& d = X_nei[facei] - X_[bCellID];

                    Ugrad[bCellID] +=
                        Ainv_[bCellID] & (U_nei[facei]-U[bCellID])*d;
                }
            }
        }

        Ugrad.correctBoundaryConditions();
    }

    tvf.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volTensorField gradientSchemes::gradient
(
    const GeometricField<vector, fvPatchField, volMesh>& U
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();

    UgradX = gradientSchemes::gradient(U.component(0));
    UgradY = gradientSchemes::gradient(U.component(1));
    UgradZ = gradientSchemes::gradient(U.component(2));

    tmp<GeometricField<tensor, fvPatchField, volMesh> > ttf
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> Ugrad = ttf();

    forAll(mesh_.cells(), cellID)
    {
        Ugrad[cellID] = tensor(UgradX[cellID], UgradY[cellID], UgradZ[cellID]);
    }

    if( Pstream::parRun() )
    {
        Ugrad.correctBoundaryConditions();
    }

    tvf.clear();
    ttf.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
    const GeometricField<tensor, fvPatchField, volMesh>& U,
    GeometricField<tensor, fvPatchField, volMesh>& UgradX,
    GeometricField<tensor, fvPatchField, volMesh>& UgradY,
    GeometricField<tensor, fvPatchField, volMesh>& UgradZ
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

    if (Pstream::parRun())
    {
        Ux.correctBoundaryConditions();
        Uy.correctBoundaryConditions();
        Uz.correctBoundaryConditions();
    }

    UgradX = gradientSchemes::gradient(Ux);
    UgradY = gradientSchemes::gradient(Uy);
    UgradZ = gradientSchemes::gradient(Uz);

    tvf.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volTensorField gradientSchemes::localGradient
(
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Unei
) const
{
    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tvft
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> Ugrad = tvft();

    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];
        const vector& dOwn = XF_[faceID] - X_[ownID];
        const vector& dNei = XF_[faceID] - X_[neiID];

        UgradX[ownID] += AinvLocal_[ownID] & ((Unei[faceID].x()-U[ownID].x())*dOwn);
        UgradY[ownID] += AinvLocal_[ownID] & ((Unei[faceID].y()-U[ownID].y())*dOwn);
        UgradZ[ownID] += AinvLocal_[ownID] & ((Unei[faceID].z()-U[ownID].z())*dOwn);

        UgradX[neiID] += AinvLocal_[neiID] & ((Unei[faceID].x()-U[neiID].x())*dNei);
        UgradY[neiID] += AinvLocal_[neiID] & ((Unei[faceID].y()-U[neiID].y())*dNei);
        UgradZ[neiID] += AinvLocal_[neiID] & ((Unei[faceID].z()-U[neiID].z())*dNei);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundary()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            vector d = XF_.boundaryField()[patchID][facei] - X_[bCellID];

            UgradX[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].x() - U[bCellID].x())*d);

            UgradY[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].y() - U[bCellID].y())*d);

            UgradZ[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].z() - U[bCellID].z())*d);

            if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            {
                const label& faceID =
                    mesh_.boundary()[patchID].start() + facei;

                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei];
                    vector d = XN_[nodeID] - X_[bCellID];

                    UgradX[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].x() - U[bCellID].x())*d);

                    UgradY[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].y() - U[bCellID].y())*d);

                    UgradZ[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].z() - U[bCellID].z())*d);

                    for (int i=0; i<7; i++)
                    {
                        d =
                            ((((i+1)*XN_[nodeID])
                          + ((7-i)*XF_.boundaryField()[patchID][facei]))/8.0)
                          - X_[bCellID];

                        UgradX[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].x() - U[bCellID].x())*d);

                        UgradY[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].y() - U[bCellID].y())*d);

                        UgradZ[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].z() - U[bCellID].z())*d);
                    }
                }
            }
        }
    }

    forAll(mesh_.cells(), cellID)
    {
        Ugrad[cellID] = tensor(UgradX[cellID], UgradY[cellID], UgradZ[cellID]);
    }

    tvf.clear();
    tvft.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvPatchField, volMesh>& Ugrad,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Um,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Up
)
{
    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];

        Um[faceID] = U[ownID] + (Ugrad[ownID] & (XF_[faceID] - X_[ownID]));
        Up[faceID]  = U[neiID] + (Ugrad[neiID] & (XF_[faceID] - X_[neiID]));
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        forAll(mesh_.boundaryMesh()[patchID],facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            U.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + ( Ugrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            Um.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + ( Ugrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Um,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Up
)
{
    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];

        Um[faceID] = U[ownID] + (Ugrad[ownID] & (XF_[faceID] - X_[ownID]));
        Up[faceID] = U[neiID] + (Ugrad[neiID] & (XF_[faceID] - X_[neiID]));
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            Um.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + (Ugrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<tensor, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& UxGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UyGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UzGrad,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Um,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Up
)
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstruct("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            U.dimensions()
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "reconstruct("+Um.name()+')',
                Um.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Um.dimensions()
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh> UmX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmZ = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpZ = tsf();

    op.decomposeTensor(Um, UmX, UmY, UmZ);
    op.decomposeTensor(Up, UpX, UpY, UpZ);

    gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX);
    gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY);
    gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ);

    forAll(own_, faceID)
    {
        Um[faceID] = tensor(UmX[faceID], UmY[faceID], UmZ[faceID]);
        Up[faceID] = tensor(UpX[faceID], UpY[faceID], UpZ[faceID]);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            const vector& reconsX =
                Ux[bCellID] + (UxGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            const vector& reconsY =
                Uy[bCellID] + (UyGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            const vector& reconsZ =
                Uz[bCellID] + (UzGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            U.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);

            Um.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);
        }
    }

    tvf.clear();
    tsf.clear();
}


//-------------------limiters-------------------------------------

inline scalar barthJespersenLimiter
(
    const scalar deltaU,        // u_{eβ} - u_e  = grad · dx
    const scalar deltaUmax,     // u_e^max - u_e
    const scalar deltaUmin      // u_e^min - u_e
)
{
    if (mag(deltaU) < SMALL)
    {
        return 1.0;
    }

    if (deltaU > 0.0)
    {
        return min(1.0, deltaUmax / deltaU);
    }
    else
    {
        return min(1.0, deltaUmin / deltaU);
    }
}

void gradientSchemes::reconstruct
(
    GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvPatchField, volMesh>& Ugrad,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Um,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Up,
    GeometricField<scalar, fvPatchField, volMesh>& phi
)
{
    // ============================================================
    // First-order fallback
    // ============================================================
    if (order_ != "second")
    {
        forAll(mesh_.cells(), cellI)
        {
            phi[cellI] = 0.0;
        }
    }
    else if (limiter_ == "no")
    {
        forAll(mesh_.cells(), cellI)
        {
            phi[cellI] = 1.0;
        }
    }
    else if (limiter_ == "yes")
    {
        // ========================================================
        // Step 1: Compute neighbourhood extrema (Algorithm 3.1-1)
        // ========================================================

        volScalarField Umin
        (
            IOobject("Umin", U.instance(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
            U
        );

        volScalarField Umax
        (
            IOobject("Umax", U.instance(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
            U
        );

        const labelUList& owner = mesh_.owner();
        const labelUList& neighbour = mesh_.neighbour();

        forAll(owner, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            Umin[own] = min(Umin[own], U[nei]);
            Umin[nei] = min(Umin[nei], U[own]);

            Umax[own] = max(Umax[own], U[nei]);
            Umax[nei] = max(Umax[nei], U[own]);
        }

        // ========================================================
        // Step 2 & 3: Compute cell limiter (Algorithm 3.1-2,3,4)
        // ========================================================

        forAll(mesh_.cells(), cellI)
        {
            scalar limiterCell = 1.0;

            const labelList& faces = mesh_.cells()[cellI];

            forAll(faces, fI)
            {
                const label faceI = faces[fI];

                vector dX = XF_[faceI] - X_[cellI];

                // Unlimited extrapolation:
                // u_{eβ} - u_e = grad · dx
                scalar deltaU =
                    (Ugrad[cellI] & dX);

                scalar deltaUmax = Umax[cellI] - U[cellI];
                scalar deltaUmin = Umin[cellI] - U[cellI];

                scalar phiFace =
                    barthJespersenLimiter(deltaU, deltaUmax, deltaUmin);

                limiterCell = min(limiterCell, phiFace);
            }

            phi[cellI] = limiterCell;
        }
    }

    // ============================================================
    // Step 4: Final MUSCL reconstruction using cell limiter
    // ============================================================

    forAll(mesh_.owner(), faceI)
    {
        const labelUList& owner = mesh_.owner();
        const labelUList& neighbour = mesh_.neighbour();

        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        Um[faceI] =
            U[own]
          + phi[own] * (Ugrad[own] & (XF_[faceI] - X_[own]));

        Up[faceI] =
            U[nei]
          + phi[nei] * (Ugrad[nei] & (XF_[faceI] - X_[nei]));
    }

    // ============================================================
    // Boundary faces
    // ============================================================

    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundaryMesh()[patchI], faceI)
        {
            const label cellI =
                mesh_.boundaryMesh()[patchI].faceCells()[faceI];

            scalar Ub =
                U[cellI]
              + phi[cellI]
              * (Ugrad[cellI]
              & (XF_.boundaryField()[patchI][faceI] - X_[cellI]));

            Um.boundaryFieldRef()[patchI][faceI] = Ub;
            U.boundaryFieldRef()[patchI][faceI]  = Ub;
        }
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void gradientSchemes::reconstruct
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Um,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Up,
    GeometricField<vector, fvPatchField, volMesh>& phi
)
{
  tmp<GeometricField<scalar, fvPatchField, volMesh> > tvf
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstruct("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            U.dimensions()
        )
    );
    GeometricField<scalar, fvPatchField, volMesh> Ux = tvf();
    GeometricField<scalar, fvPatchField, volMesh> Uy = tvf();
    GeometricField<scalar, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    Ux= U.component(0);
    Uy= U.component(1);
    Uz= U.component(2);

      tmp<GeometricField<scalar, fvPatchField, volMesh> > tvfPhi
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstructPhi("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            phi.dimensions()
        )
    );

    GeometricField<scalar, fvPatchField, volMesh> phiX = tvfPhi();
    GeometricField<scalar, fvPatchField, volMesh> phiY = tvfPhi();
    GeometricField<scalar, fvPatchField, volMesh> phiZ = tvfPhi();

    phiX= phi.component(0);
    phiY= phi.component(1);
    phiZ= phi.component(2);




  tmp<GeometricField<vector, fvPatchField, volMesh> > tvf1
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstruct("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Ugrad.dimensions()
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UxGrad = tvf1();
    GeometricField<vector, fvPatchField, volMesh> UyGrad = tvf1();
    GeometricField<vector, fvPatchField, volMesh> UzGrad = tvf1();

    op.decomposeTensor(Ugrad, UxGrad, UyGrad, UzGrad);



    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<scalar, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "reconstruct("+Um.name()+')',
                Um.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Um.dimensions()
        )
    );
    GeometricField<scalar, fvsPatchField, surfaceMesh> UmX = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UmY = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UmZ = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UpX = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UpY = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UpZ = tsf();



    // Info << "------vector Reconstruct---------"<<endl;

    gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX, phiX);
    gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY, phiY);
    gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ, phiZ);
    // gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX);
    // gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY);
    // gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ);

    forAll(own_, faceID)
    {
        Um[faceID] = vector(UmX[faceID], UmY[faceID], UmZ[faceID]);
        Up[faceID] = vector(UpX[faceID], UpY[faceID], UpZ[faceID]);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            const scalar& reconsX =
                Ux[bCellID] + ((UxGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID])));

            const scalar& reconsY =
                Uy[bCellID] + ((UyGrad[bCellID]
               & (XF_.boundaryField()[patchID][facei] - X_[bCellID])));

            const scalar& reconsZ =
                Uz[bCellID] + ( (UzGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID])));

            U.boundaryFieldRef()[patchID][facei] =
                vector(reconsX, reconsY, reconsZ);

            Um.boundaryFieldRef()[patchID][facei] =
                vector(reconsX, reconsY, reconsZ);
        }
    }

    tvf.clear();
    tsf.clear();
    tvf1.clear();
    tvfPhi.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<tensor, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& UxGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UyGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UzGrad,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Um,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Up,
    GeometricField<tensor, fvPatchField, volMesh>& phi

)
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstruct("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            U.dimensions()
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

        tmp<GeometricField<vector, fvPatchField, volMesh> > tvfphi
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstruct("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            phi.dimensions()
        )
    );

    GeometricField<vector, fvPatchField, volMesh> phiX = tvfphi();
    GeometricField<vector, fvPatchField, volMesh> phiY = tvfphi();
    GeometricField<vector, fvPatchField, volMesh> phiZ = tvfphi();
    op.decomposeTensor(phi, phiX, phiY, phiZ);

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "reconstruct("+Um.name()+')',
                Um.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Um.dimensions()
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh> UmX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmZ = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpZ = tsf();

    op.decomposeTensor(Um, UmX, UmY, UmZ);
    op.decomposeTensor(Up, UpX, UpY, UpZ);

    // Info << "------Tendor Reconstruct---------"<<endl;
    gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX, phiX);
    gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY, phiY);
    gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ, phiZ);
    // gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX);
    // gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY);
    // gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ);

    forAll(mesh_.cells(), cellID)
    {
        phi[cellID] = tensor(phiX[cellID], phiY[cellID], phiZ[cellID]);
    }

    forAll(own_, faceID)
    {
        Um[faceID] = tensor(UmX[faceID], UmY[faceID], UmZ[faceID]);
        Up[faceID] = tensor(UpX[faceID], UpY[faceID], UpZ[faceID]);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];


            const vector& reconsX =
                Ux[bCellID] + (  (UxGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID])));

            const vector& reconsY =
                Uy[bCellID] + ( (UyGrad[bCellID]
               & (XF_.boundaryField()[patchID][facei] - X_[bCellID])));

            const vector& reconsZ =
                Uz[bCellID] + (  (UzGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID])));

            U.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);

            Um.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);
        }
    }

    tvf.clear();
    tsf.clear();
    tvfphi.clear();
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
