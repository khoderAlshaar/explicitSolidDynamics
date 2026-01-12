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

Class
    Foam::sinusoidalTractionLinearMomentumFvPatchVectorField

Description
    Specifies contact linear momentum for a traction boundary patch according
    to the Riemann solver. External boundary traction or pressure can also
    be specified.

Usage
    \table
        Property     | Description             | Required    | Default value
        loadingType  | nature of loading: none,| no          | none
                     | traction or pressure_   |             |
        traction     | traction vector         | no          | vector(0,0,0)
        pressure     | pressure                | no          | 0.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    free
    {
        type           tractionLinearMomentum;
        loadingType    traction;
        traction       (0 10 0);
    }
    free
    {
        type           tractionLinearMomentum;
        loadingType    pressure;
        pressure       5.0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    sinusoidalTractionLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef sinusoidalTractionLinearMomentumFvPatchVectorField_H
#define sinusoidalTractionLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class sinusoidalTractionLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sinusoidalTractionLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

private:

    // Sinusoidal traction parameters
    scalar K_;
    scalar omega_;
    scalar phase_;
    vector direction_;


public:

    // Runtime type information
    TypeName("sinusoidalTractionLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        sinusoidalTractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sinusoidalTractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given sinusoidalTractionLinearMomentumFvPatchVectorField
        //  onto a new patch
        sinusoidalTractionLinearMomentumFvPatchVectorField
        (
            const sinusoidalTractionLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sinusoidalTractionLinearMomentumFvPatchVectorField
        (
            const sinusoidalTractionLinearMomentumFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new sinusoidalTractionLinearMomentumFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sinusoidalTractionLinearMomentumFvPatchVectorField
        (
            const sinusoidalTractionLinearMomentumFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new sinusoidalTractionLinearMomentumFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            // virtual word& loading()
            // {
            //     return loadingType_;
            // }

            // virtual vector& traction()
            // {
            //     return t_P_;
            // }

            // virtual scalar& pressure()
            // {
            //     return p_P_;
            // }

        // Mapping functions

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

        // Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //