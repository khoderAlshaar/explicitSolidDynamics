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

Class
    analyticalPlateHoleTractionTractionFvPatchVectorField

Description
    Traction boundary condition for plateHole test case to enforce
    analytical traction at the boundaries

Author
   Philip Cardiff, UCD.
   Zeljko Tukovic, Zagreb.

SourceFiles
    analyticalPlateHoleTractionTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalPlateHoleTractionTractionFvPatchVectorField_H
#define analyticalPlateHoleTractionTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class analyticalPlateHoleTractionTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class analyticalPlateHoleTractionTractionFvPatchVectorField
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

    //- Runtime type information
    TypeName("analyticalPlateHoleTractionTraction");


    // Constructors

        //- Construct from patch and internal field
        analyticalPlateHoleTractionTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        analyticalPlateHoleTractionTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  analyticalPlateHoleTractionTractionFvPatchVectorField onto a new patch
        analyticalPlateHoleTractionTractionFvPatchVectorField
        (
            const analyticalPlateHoleTractionTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        analyticalPlateHoleTractionTractionFvPatchVectorField
        (
            const analyticalPlateHoleTractionTractionFvPatchVectorField&
// #ifdef OPENFOAM_ORG
        ) = delete;
// #else
//         );

//         //- Construct and return a clone
//         virtual tmp<fvPatchVectorField> clone() const
//         {
//             return tmp<fvPatchVectorField>
//             (
//                 new analyticalPlateHoleTractionTractionFvPatchVectorField(*this)
//             );
//         }
// #endif

        //- Construct as copy setting internal field reference
        analyticalPlateHoleTractionTractionFvPatchVectorField
        (
            const analyticalPlateHoleTractionTractionFvPatchVectorField&,
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
                new analyticalPlateHoleTractionTractionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
