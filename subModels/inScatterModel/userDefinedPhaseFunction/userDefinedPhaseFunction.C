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

#include "userDefinedPhaseFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(userDefinedPhaseFunction, 0);

        addToRunTimeSelectionTable
        (
            inScatterModel,
            userDefinedPhaseFunction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::userDefinedPhaseFunction::userDefinedPhaseFunction
(
    const dictionary& dict
)
:
    inScatterModel(dict),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs"))
{
    const dictionary& functionDicts = dict.optionalSubDict(typeName + "Coeffs");
    functionDicts.lookup("inScatter") >> inScatter_;
    functionDicts.lookup("nBand") >> nBands_;
    functionDicts.lookup("phaseFunctionAngleNum") >> phaseFunctionAngleNum_;
    phaseFunc_.setSize(nBands_*phaseFunctionAngleNum_);   
    deltaAngle = Foam::constant::mathematical::pi/scalar(phaseFunctionAngleNum_-1 );   
     
	if(inScatter_)
    {
		functionDicts.lookup("phaseFunction") >> phaseFunc_;
	//	Info << "phaseFunc_  " << phaseFunc_ << endl;
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::light::userDefinedPhaseFunction::correct
(
    const scalar rayCos,
    const label  iBand
)const
{
    scalar angle = Foam::acos(rayCos);
    label i0 = label (angle / deltaAngle + 0.2);

    if (i0 >= 0 && i0 < phaseFunctionAngleNum_)
    {
        return phaseFunc_[i0 + iBand * phaseFunctionAngleNum_];
    }
    else
    {
        return 0.0;
    }
}


// ************************************************************************* //
