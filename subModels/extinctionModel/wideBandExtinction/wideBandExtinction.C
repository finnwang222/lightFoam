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

\*---------------------------------------------------------------------------*/

#include "wideBandExtinction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(wideBandExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            wideBandExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::wideBandExtinction::wideBandExtinction
(
    const dictionary& dict
)
:
    extinctionModel(dict),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    totalWaveLength_(0)
{
    const dictionary& functionDicts = dict.optionalSubDict(typeName + "Coeffs");
    functionDicts.lookup("absorption") >> absorption_;
    functionDicts.lookup("scattering") >> scattering_;
    functionDicts.lookup("totalWaveLength") >> totalWaveLength_;
    functionDicts.lookup("nBand") >> nBands_;
    functionDicts.lookup("bioDensity") >> bioDensity_;
    aBand_.setSize(nBands_);
    sBand_.setSize(nBands_);
    
    if(absorption_)
    {
        functionDicts.lookup("absorptionCoeff") >> aBand_;
        forAll (aBand_, i) aBand_[i] = aBand_[i] * bioDensity_;
    }
    else
    {
        forAll (aBand_, i) aBand_[i] = 0.0;
    }

    if(scattering_)
    {
        functionDicts.lookup("scatteringCoeff") >> sBand_;
        forAll (sBand_, i) sBand_[i] = sBand_[i] * bioDensity_;
    }
    else
    {
        forAll (sBand_, i) sBand_[i] = 0.0;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::wideBandExtinction::~wideBandExtinction()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::light::wideBandExtinction::a(const label bandI) const
{
    return aBand_[bandI];
}


Foam::scalar
Foam::light::wideBandExtinction::s(const label bandI) const
{
    return sBand_[bandI];
}

// ************************************************************************* //
