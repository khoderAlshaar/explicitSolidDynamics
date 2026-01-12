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
#include "sinusoidalTractionLinearMomentumFvPatchVectorField.H"
// #include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sinusoidalTractionLinearMomentumFvPatchVectorField::
sinusoidalTractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    K_(1e-3),
    omega_(constant::mathematical::pi/20.0),
    phase_(-constant::mathematical::pi/2.0),
    direction_(vector(1,0,0))
{}


sinusoidalTractionLinearMomentumFvPatchVectorField::
sinusoidalTractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    K_(dict.lookupOrDefault("K", 1e-3)),
    omega_(dict.lookupOrDefault("omega", constant::mathematical::pi/20.0)),
    phase_(dict.lookupOrDefault("phase", -constant::mathematical::pi/2.0)),
    direction_(dict.lookupOrDefault("direction", vector(1,0,0)))
{
    direction_ /= mag(direction_);   // normalize
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


sinusoidalTractionLinearMomentumFvPatchVectorField::
sinusoidalTractionLinearMomentumFvPatchVectorField
(
    const sinusoidalTractionLinearMomentumFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    K_(ptf.K_),
    omega_(ptf.omega_),
    phase_(ptf.phase_),
    direction_(ptf.direction_)
{}


sinusoidalTractionLinearMomentumFvPatchVectorField::
sinusoidalTractionLinearMomentumFvPatchVectorField
(
    const sinusoidalTractionLinearMomentumFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    K_(ptf.K_),
    omega_(ptf.omega_),
    phase_(ptf.phase_),
    direction_(ptf.direction_)
{}


sinusoidalTractionLinearMomentumFvPatchVectorField::
sinusoidalTractionLinearMomentumFvPatchVectorField
(
    const sinusoidalTractionLinearMomentumFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    K_(ptf.K_),
    omega_(ptf.omega_),
    phase_(ptf.phase_),
    direction_(ptf.direction_)
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sinusoidalTractionLinearMomentumFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    sinusoidalTractionLinearMomentumFvPatchVectorField::rmap(ptf, addr);
}


void sinusoidalTractionLinearMomentumFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M =
        patch().lookupPatchField<surfaceVectorField, vector>("lm_M");

    const fvsPatchField<vector>& t_M =
        patch().lookupPatchField<surfaceVectorField, vector>("t_M");

    const fvsPatchField<tensor>& S_t =
        patch().lookupPatchField<surfaceTensorField, tensor>("S_t");

    const scalar t = this->db().time().value();

    scalar amplitude =
        K_*(sin(omega_*t + phase_) + 1.0);

    vector traction = amplitude * direction_;

    fvsPatchField<vector> lm_C(lm_M);

     lm_C = lm_M + (S_t & (traction - t_M));

    operator==(lm_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void sinusoidalTractionLinearMomentumFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("K") << K_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << omega_ << token::END_STATEMENT << nl;
    os.writeKeyword("phase") << phase_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction") << direction_ << token::END_STATEMENT << nl;
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchVectorField,
    sinusoidalTractionLinearMomentumFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //