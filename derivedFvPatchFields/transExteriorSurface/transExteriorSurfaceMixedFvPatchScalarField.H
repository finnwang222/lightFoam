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
    transExteriorSurfaceMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature,
    to be used by the conjugate heat trnasfer solver. 
    Both sides use a mix of zeroGradient and neigbourvalue.

SourceFiles
    transExteriorSurfaceMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef transExteriorSurfaceMixedFvPatchScalarField_H
#define transExteriorSurfaceMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

/*---------------------------------------------------------------------------*\
    Class transExteriorSurfaceMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class transExteriorSurfaceMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private Data

        //- Name of field on the neighbour region
        
        scalar nNbg_;  //index of refraction
        scalar nOwn_;
        scalar I0_;    //incident light intensity
        label  nBands_;
        scalar diffuseFraction_;
        scalar beamWidthPhi_;
        scalar beamWidthTheta_;
        vector beamDir_;  //incident light direction
        bool beamNormToSurf_;
        label initFlag_;
        DynamicList<scalar> lightBandDist_;
        autoPtr<scalarField> spectacularRefraction;
        autoPtr<scalarField> diffusiveRefraction;
        scalar beamAnglePhi_;
        scalar beamAngleTheta_;

public:

    //- Runtime type information
    TypeName("transExteriorSurface");


    // Constructors

        //- Construct from patch and internal field
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  transExteriorSurfaceMixedFvPatchScalarField onto a
        //  new patch
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const transExteriorSurfaceMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        transExteriorSurfaceMixedFvPatchScalarField
        (
            const transExteriorSurfaceMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new transExteriorSurfaceMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const transExteriorSurfaceMixedFvPatchScalarField&,
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
                new transExteriorSurfaceMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        virtual void dirToAngle(const vector&, scalar&, scalar&) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
