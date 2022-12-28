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
    Foam::light::userDefinedPhaseFunction

Description
    User define different phase function for different media

SourceFiles
    inScatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef userDefinedPhaseFunction_H
#define userDefinedPhaseFunction_H
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
                       Class userDefinedPhaseFunction Declaration
\*---------------------------------------------------------------------------*/

class userDefinedPhaseFunction
:
    public inScatterModel
{
public:

    DynamicList<scalar> phaseFunc_;

    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        Switch inScatter_;
        label  nBands_;
        label  phaseFunctionAngleNum_;
        scalar deltaAngle;



public:

    //- Runtime type information
    TypeName("userDefinedPhaseFunction");


    // Constructors

        //- Construct from components
        userDefinedPhaseFunction
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~userDefinedPhaseFunction() = default;


    // Member Functions

        inline bool inScatter() const
        {
            return inScatter_;
        }

        Foam::scalar correct
        (
            const scalar rayCos,
            const label  iBand
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
