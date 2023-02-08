/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    multiBandTransInteriorCoupledFvPatchScalarField

Description
    Mixed boundary condition for temperature,
    to be used by the conjugate heat trnasfer solver. 
    Both sides use a mix of zeroGradient and neigbourvalue.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                    multiBandTransInteriorCoupled;
        neighbourFieldName      T;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

SourceFiles
    multiBandTransInteriorCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef multiBandTransInteriorCoupledFvPatchScalarField_H
#define multiBandTransInteriorCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

/*---------------------------------------------------------------------------*\
    Class multiBandTransInteriorCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class multiBandTransInteriorCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private Data

        //- Name of field on the neighbour region
        
        scalar nNbg_;  //index of refraction
        scalar nOwn_;
        scalar diffuseFraction_;
        label nBands_;
        DynamicList<scalar> lightBandDist_;

public:

    //- Runtime type information
    TypeName("multiBandTransInteriorCoupled");


    // Constructors

        //- Construct from patch and internal field
        multiBandTransInteriorCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        multiBandTransInteriorCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  multiBandTransInteriorCoupledFvPatchScalarField onto a
        //  new patch
        multiBandTransInteriorCoupledFvPatchScalarField
        (
            const multiBandTransInteriorCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new multiBandTransInteriorCoupledFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        multiBandTransInteriorCoupledFvPatchScalarField
        (
            const multiBandTransInteriorCoupledFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new multiBandTransInteriorCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
