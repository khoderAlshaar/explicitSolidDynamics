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

#include "analyticalPlateHoleTractionLinearMomentumFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
// #include "mechanicalModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "coordinateSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

symmTensor analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::plateHoleSolution
(
    const vector& C
)
{
    tensor sigma = tensor::zero;

    // Calculate radial coordinate
    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));

    // Calculate circumferential coordinate
    scalar theta = Foam::atan2(C.y(), C.x());

    coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

    sigma.xx() =
        T_*(1 - sqr(holeR_)/sqr(r))/2
      + T_
       *(1 + 3*pow(holeR_,4)/pow(r,4) - 4*sqr(holeR_)/sqr(r))*::cos(2*theta)/2;

    sigma.xy() =
      - T_
       *(1 - 3*pow(holeR_,4)/pow(r,4) + 2*sqr(holeR_)/sqr(r))*::sin(2*theta)/2;

    sigma.yx() = sigma.xy();

    sigma.yy() =
        T_*(1 + sqr(holeR_)/sqr(r))/2
      - T_*(1 + 3*pow(holeR_,4)/pow(r,4))*::cos(2*theta)/2;


    // Transformation to Cartesian coordinate system
// #ifdef OPENFOAM_ORG
    sigma = ((cs.R().R() & sigma) & cs.R().R().T());
// #else
//     sigma = ((cs.R() & sigma) & cs.R().T());
// #endif

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::
analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    T_(0.0),
    holeR_(0.0)
{}


analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::
analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius")))
{
    updateCoeffs();
}


analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::
analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
(
    const analyticalPlateHoleTractionLinearMomentumFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    T_(ptf.T_),
    holeR_(ptf.holeR_)
{}


analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::
analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
(
    const analyticalPlateHoleTractionLinearMomentumFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    T_(rifvpvf.T_),
    holeR_(rifvpvf.holeR_)
{}


analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::
analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
(
    const analyticalPlateHoleTractionLinearMomentumFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    T_(rifvpvf.T_),
    holeR_(rifvpvf.holeR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("lm_M");

    const fvsPatchField<vector>& t_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("t_M");

    const fvsPatchField<tensor>& S_t_ =
        patch().lookupPatchField<surfaceTensorField, tensor>("S_t");

    fvsPatchField<vector> lm_C(lm_M_);

    fvsPatchField<vector> t_C(lm_M_);

    vectorField& trac = t_C;
        // Patch unit normals
    vectorField n(patch().nf());

    // Patch face centres
    const vectorField& Cf = patch().Cf();

    // const volTensorField& sigma =
    //             mesh.lookupObject<volTensorField>("P");

    forAll(trac, faceI)
    {
        vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);
        vector curN = n[faceI];

        if (patch().name() == "hole")
        {
            curC /= mag(curC);
            curC *= holeR_;

            curN = -curC/mag(curC);
        }

        trac[faceI] = (n[faceI] & plateHoleSolution(curC));
    }


        lm_C = lm_M_ + (S_t_ & ((trac) - t_M_));



    this->operator==(lm_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void analyticalPlateHoleTractionLinearMomentumFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    writeEntry(os, "value", *this);

        os.writeKeyword("farFieldTractionX")
        << T_ << token::END_STATEMENT << nl;

    os.writeKeyword("holeRadius")
        << holeR_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //