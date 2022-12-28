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
    Foam::light::wideBandExtinction


Description


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef lightWideBandExtinction_H
#define lightWideBandExtinction_H

#include "extinctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

/*---------------------------------------------------------------------------*\
                  Class wideBandExtinction Declaration
\*---------------------------------------------------------------------------*/

class wideBandExtinction
:
    public extinctionModel
{
    // Private data

        //- Extinction model dictionary
        dictionary coeffsDict_;

        Switch absorption_;

        Switch scattering_;

        scalar totalWaveLength_;

        label nBands_;

        scalar bioDensity_;

        DynamicList<scalar>   aBand_;

        DynamicList<scalar>   sBand_;


public:

    //- Runtime type information
    TypeName("wideBandExtinction");


    // Constructors

        //- Construct from components
        wideBandExtinction(const dictionary& dict);


    //- Destructor
    virtual ~wideBandExtinction();


    // Member Functions

        scalar a(const label bandI) const;

        scalar s(const label bandI) const;

        inline scalar totalWaveLength() const
        {
            return totalWaveLength_;
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
