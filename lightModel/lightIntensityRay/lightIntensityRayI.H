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

inline const Foam::volScalarField&
Foam::light::lightIntensityRay::I() const
{
    return I_;
}

inline const Foam::vector& Foam::light::lightIntensityRay::d() const
{
    return d_;
}


inline const Foam::vector& Foam::light::lightIntensityRay::dAve() const
{
    return dAve_;
}

inline Foam::label Foam::light::lightIntensityRay::iBand() const
{
    return iBand_;
}

inline Foam::label Foam::light::lightIntensityRay::iAngle() const
{
    return iAngle_;
}


inline Foam::scalar Foam::light::lightIntensityRay::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::light::lightIntensityRay::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::light::lightIntensityRay::omega() const
{
    return omega_;
}

// ************************************************************************* //
