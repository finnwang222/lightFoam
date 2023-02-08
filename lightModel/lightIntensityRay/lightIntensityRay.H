/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::light::lightIntensityRay

Description
    Light intensity for a ray in a given direction

SourceFiles
    lightIntensityRay.C

\*---------------------------------------------------------------------------*/

#ifndef lightIntensityRay_H
#define lightIntensityRay_H

#include "extinctionModel.H"
#include "inScatterModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

// Forward declaration of classes
class lightDOM;

/*---------------------------------------------------------------------------*\
                    Class lightIntensityRay Declaration
\*---------------------------------------------------------------------------*/

class lightIntensityRay
{
public:

    static const word intensityPrefix;


private:

    // Private data

        //- Reference to the owner lightDOM object
        const lightDOM& dom_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Total radiative intensity / [W/m2]
        autoPtr<volScalarField> I_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Number of wavelengths/bands
        label iBand_;

        label iAngle_;
        
        //- Global ray id - incremented in constructor
        static label rayId;

    // Private Member Functions

        //- No copy construct
        lightIntensityRay(const lightIntensityRay&) = delete;

        //- No copy assignment
        void operator=(const lightIntensityRay&) = delete;


public:

    // Constructors

        //- Construct form components
        lightIntensityRay
        (
            const lightDOM& dom,
            const fvMesh& mesh,
            const label iBand,
            const label iAngle,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta
        );


    //- Destructor
    ~lightIntensityRay();


    // Member functions

        // Edit

            //- Update radiative intensity on i direction
            scalar correct();

            //- Initialise the ray in i direction
            void init
            (
                const scalar phi,
                const scalar theta,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

        // Access

            //- Return intensity
            inline const volScalarField& I() const;

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            //- Return the number of bands
            inline label iBand() const;

            inline label iAngle() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            void updateBoundary();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lightIntensityRayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
