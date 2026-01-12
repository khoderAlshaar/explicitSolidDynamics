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
    Foam::analyticalPlateHoleTractionLinearMomentumFvPatchVectorField

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
    analyticalPlateHoleTractionLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalPlateHoleTractionLinearMomentumFvPatchVectorField_H
#define analyticalPlateHoleTractionLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class analyticalPlateHoleTractionLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

private:

    // Private Data

        //- Far field traction
        //  It is assumed that it acts in the positive X direction
        const scalar T_;

        //- Hole radius
        const scalar holeR_;


    // Private Member Functions

        //- Analytical solution
        symmTensor plateHoleSolution(const vector& C);


public:

    // Runtime type information
    TypeName("analyticalPlateHoleTractionLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
        //  onto a new patch
        analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
        (
            const analyticalPlateHoleTractionLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
        (
            const analyticalPlateHoleTractionLinearMomentumFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new analyticalPlateHoleTractionLinearMomentumFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        analyticalPlateHoleTractionLinearMomentumFvPatchVectorField
        (
            const analyticalPlateHoleTractionLinearMomentumFvPatchVectorField&,
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
                new analyticalPlateHoleTractionLinearMomentumFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

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