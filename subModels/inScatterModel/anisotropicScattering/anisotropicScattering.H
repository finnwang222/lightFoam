/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::light::anisotropicScattering

Description
    Anisotropic light scatter coefficient

SourceFiles
    inScatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropicScattering_H
#define anisotropicScattering_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "inScatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

/*---------------------------------------------------------------------------*\
                       Class anisotropicScattering Declaration
\*---------------------------------------------------------------------------*/

class anisotropicScattering
:
    public inScatterModel
{
public:

    DynamicList<scalar> phaseFuncCoef_;

    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        label nBands_;


public:

    //- Runtime type information
    TypeName("anisotropicScattering");


    // Constructors

        //- Construct from components
        anisotropicScattering
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~anisotropicScattering();


    // Member Functions

        Foam::scalar correct
        (
            const scalar rayCos,
            const label  iBand
        )const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
