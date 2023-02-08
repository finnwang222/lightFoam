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
    Foam::radiation::absorptionEmissionModel

Description
    Model to supply absorption and emission coefficients for radiation
    modelling

\*---------------------------------------------------------------------------*/

#ifndef lightExtinctionModel_H
#define lightExtinctionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

/*---------------------------------------------------------------------------*\
                   Class extinctionModel Declaration
\*---------------------------------------------------------------------------*/

class extinctionModel
{

protected:

    // Protected data

        //- Radiation model dictionary
        const dictionary dict_;


public:

    //- Runtime type information
    TypeName("extinctionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            extinctionModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from components
        extinctionModel
        (
            const dictionary& dict
        );


    //- Selector
    static autoPtr<extinctionModel> New
    (
        const dictionary& dict
    );


    //- Destructor
    virtual ~extinctionModel();


    // Member Functions

        // Access

            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }


            // Extinction coefficient


                virtual scalar a(const label bandI = 0) const;

                virtual scalar k(const label bandI = 0) const;

                virtual scalar s(const label bandI = 0) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
