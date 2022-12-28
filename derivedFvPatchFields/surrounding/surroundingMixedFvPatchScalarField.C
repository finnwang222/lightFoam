/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "surroundingMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "lightDOM.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    surroundLight_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const surroundingMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    surroundLight_(ptf.surroundLight_)
{}


surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    surroundLight_(readScalar(dict.lookup("surroundLight")))
{

    if (dict.found("refValue"))
    {
         fvPatchScalarField::operator=
         (
            scalarField("value", dict, p.size())
         );

        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
        
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchField::operator = (refValue());
    }
}


surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const surroundingMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    surroundLight_(ptf.surroundLight_)
{}


surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const surroundingMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf),
    surroundLight_(ptf.surroundLight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void surroundingMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    scalarField& Iw = *this;

    const lightModel& light = db().lookupObject<lightModel>("lightProperties");
    const lightDOM& dom(refCast<const lightDOM>(light));

    label rayId = -1;
    dom.setRayId(internalField().name(), rayId);

    const vectorField& n = patch().Sf()/patch().magSf();
    const vector& bdRayDir = dom.IRay(rayId).d();
    const scalar bdOmega = dom.IRay(rayId).omega();
    const label nBand = dom.nBand();

    forAll(Iw, faceI)
    {
        if((-n[faceI] & bdRayDir) > 0.0)
        {
            refValue()[faceI] = surroundLight_/nBand/4/pi*bdOmega;
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
        }
        else
        {
            refValue()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 0.0;
        }
    }

    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs();
}


void surroundingMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("surroundingLight ") << surroundLight_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    surroundingMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam


// ************************************************************************* //
