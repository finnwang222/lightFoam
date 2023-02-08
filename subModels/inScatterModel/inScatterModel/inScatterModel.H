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
    Foam::light::inScatterModel

Description
    Base class for light scattering

\*---------------------------------------------------------------------------*/

#ifndef inScatterModel_H
#define inScatterModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

/*---------------------------------------------------------------------------*\
                        Class scatterModel Declaration
\*---------------------------------------------------------------------------*/

class inScatterModel
{

protected:

    // Protected data

        const dictionary dict_;

public:

    //- Runtime type information
    TypeName("inScatterModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            inScatterModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from components
        inScatterModel
        (
            const dictionary& dict
        );


    // Selector

        static autoPtr<inScatterModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~inScatterModel();


    // Member Functions

        inline const dictionary& dict() const
        {
            return dict_;
        }

        virtual bool inScatter() const;

        virtual scalar correct
        (
            scalar angle,
            label  iBand
        )const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
