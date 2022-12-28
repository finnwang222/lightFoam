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
    Foam::light::deltaEddingtonModel

Description
    Delta Eddington Model

SourceFiles
    inScatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef deltaEddingtonModel_H
#define deltaEddingtonModel_H
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
                       Class deltaEddingtonModel Declaration
\*---------------------------------------------------------------------------*/

class deltaEddingtonModel
:
    public inScatterModel
{
public:

    DynamicList<scalar> phaseFuncCoef_;
    DynamicList<scalar> forwardScatterFactor_;

    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        label nBands_;


public:

    //- Runtime type information
    TypeName("deltaEddingtonModel");


    // Constructors

        //- Construct from components
        deltaEddingtonModel
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~deltaEddingtonModel();


    // Member Functions

        Foam::scalar correct
        (
            const scalar rayCos,
            const label iBand
        )const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
