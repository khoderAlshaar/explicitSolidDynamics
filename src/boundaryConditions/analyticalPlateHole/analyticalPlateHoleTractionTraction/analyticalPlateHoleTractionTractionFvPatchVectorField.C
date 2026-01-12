/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "analyticalPlateHoleTractionTractionFvPatchVectorField.H"
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

symmTensor analyticalPlateHoleTractionTractionFvPatchVectorField::plateHoleSolution
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

analyticalPlateHoleTractionTractionFvPatchVectorField::
analyticalPlateHoleTractionTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    T_(0.0),
    holeR_(0.0)
{}


analyticalPlateHoleTractionTractionFvPatchVectorField::
analyticalPlateHoleTractionTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius")))
{}


analyticalPlateHoleTractionTractionFvPatchVectorField::
analyticalPlateHoleTractionTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(stpvf, p, iF, mapper),
    T_(stpvf.T_),
    holeR_(stpvf.holeR_)
{}

// #ifndef OPENFOAM_ORG
// analyticalPlateHoleTractionTractionFvPatchVectorField::
// analyticalPlateHoleTractionTractionFvPatchVectorField
// (
//     const analyticalPlateHoleTractionTractionFvPatchVectorField& stpvf
// )
// :
//     fixedValueFvPatchVectorField(stpvf),
//     T_(stpvf.T_),
//     holeR_(stpvf.holeR_)
// {}
// #endif

analyticalPlateHoleTractionTractionFvPatchVectorField::
analyticalPlateHoleTractionTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(stpvf, iF),
    T_(stpvf.T_),
    holeR_(stpvf.holeR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void analyticalPlateHoleTractionTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void analyticalPlateHoleTractionTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void analyticalPlateHoleTractionTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Patch unit normals
    vectorField n(patch().nf());

    // Patch face centres
    const vectorField& Cf = patch().Cf();

    // Set the patch traction
    const fvsPatchField<vector>& lm_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("lm_M");

    fvsPatchField<vector> t_C(lm_M_);

    vectorField& trac = t_C;

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
    
    this->operator==(trac);
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void analyticalPlateHoleTractionTractionFvPatchVectorField::write(Ostream& os) const
{
    // fixedValueFvPatchVectorField::write(os);



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
    analyticalPlateHoleTractionTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
