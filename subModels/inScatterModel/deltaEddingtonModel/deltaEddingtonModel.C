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

#include "deltaEddingtonModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(deltaEddingtonModel, 0);

        addToRunTimeSelectionTable
        (
            inScatterModel,
            deltaEddingtonModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::deltaEddingtonModel::deltaEddingtonModel
(
    const dictionary& dict
)
:
    inScatterModel(dict),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs"))
{
    coeffsDict_.lookup("nBand") >> nBands_;
    phaseFuncCoef_.setSize(nBands_);
    coeffsDict_.lookup("phaseFuncCoef") >> phaseFuncCoef_;
    coeffsDict_.lookup("forwardScatteringFactor") >> forwardScatterFactor_;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::deltaEddingtonModel::~deltaEddingtonModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::light::deltaEddingtonModel::correct
(
    const scalar rayCos,
    const label iBand
)const
{
    scalar y;
    if(rayCos == 1)
    {
        y = 2 * forwardScatterFactor_[iBand] + 
        (1 - forwardScatterFactor_[iBand]) * (1 +  phaseFuncCoef_[iBand]);
    }
    else
    {
        y = (1 - forwardScatterFactor_[iBand]) * 
        (1 + phaseFuncCoef_[iBand]*rayCos );
    }

    return y;
}


// ************************************************************************* //
