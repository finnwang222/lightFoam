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

#include "anisotropicScattering.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(anisotropicScattering, 0);

        addToRunTimeSelectionTable
        (
            inScatterModel,
            anisotropicScattering,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::anisotropicScattering::anisotropicScattering
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
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::anisotropicScattering::~anisotropicScattering()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::light::anisotropicScattering::correct
(
    const scalar rayCos,
    const label  iBand
)const
{
    scalar y = 1 + phaseFuncCoef_[iBand] * rayCos;
    return y;
}


// ************************************************************************* //
