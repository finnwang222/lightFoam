/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "lightModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::light::lightModel>
Foam::light::lightModel::New
(
    const volScalarField& intensity
)
{
    word modelType("none");

    dictionary dict;

    IOobject io
    (
        "lightProperties",
        intensity.time().constant(),
        intensity.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false // Do not register
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary propDict(io);

        dict = std::move(propDict);

        dict.readEntry("lightModel", modelType);
    }
    else
    {
        Info<< "Light model not active: lightProperties not found"
            << endl;
    }

    Info<< "Selecting lightModel " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "lightModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<lightModel>(ctorPtr(intensity));
}









// ************************************************************************* //
