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

#include "reflectiveSurfaceMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "constants.H"
#include "lightDOM.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    diffuseFraction_(0.0),
    reflectionCoef_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const reflectiveSurfaceMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    diffuseFraction_(ptf.diffuseFraction_),
    reflectionCoef_(ptf.reflectionCoef_)
{}


reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction"))),
    reflectionCoef_(readScalar(dict.lookup("reflectionCoef")))
{
    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator = 
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


reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const reflectiveSurfaceMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    diffuseFraction_(ptf.diffuseFraction_),
    reflectionCoef_(ptf.reflectionCoef_)
{}


reflectiveSurfaceMixedFvPatchScalarField::
reflectiveSurfaceMixedFvPatchScalarField
(
    const reflectiveSurfaceMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    diffuseFraction_(ptf.diffuseFraction_),
    reflectionCoef_(ptf.reflectionCoef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reflectiveSurfaceMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;


    scalarField& Iw = *this;
    const lightModel& light = db().lookupObject<lightModel>("lightProperties");
    const lightDOM& dom(refCast<const lightDOM>(light));

    label rayId = -1;
    dom.setRayId(internalField().name(), rayId);

    // Get the coupling information from the mappedPatchBase
    const label patchI = patch().index();
    const vectorField n(patch().Sf() / patch().magSf());
    const label iBand = dom.IRay(rayId).iBand();
    const label nAngle = dom.nAngle();
    const label nTheta = dom.nTheta();
    const label nPhi = dom.nPhi();
    const vector& bdRayDir = dom.IRay(rayId).d();
    const scalar bdOmega = dom.IRay(rayId).omega();

    const scalar deltaPhi = pi / (2.0 * nPhi);
    const scalar deltaTheta = pi / nTheta;
    label npPhi = dom.NumPixelPhi();
    label npTheta = dom.NumPixelTheta();
    const scalar bdRayPhi = dom.IRay(rayId).phi();
    const scalar bdRayTheta = dom.IRay(rayId).theta();

    if (internalField().mesh().nSolutionD() == 2)
    {
        npTheta = 1;
    }
    if (internalField().mesh().nSolutionD() == 1)
    {
        npTheta = 1;
        npPhi = 1;
    }

    scalar spectacular = 0.0;
    scalar diffusive = 0.0;

    forAll (Iw, faceI)
    {
        spectacular = 0.0;
        diffusive = 0.0;
        vector surfNorm = -n[faceI];
        scalar cosA = surfNorm& bdRayDir;

        if (cosA > 0.0)
        {
            if (diffuseFraction_ > 0)
            {
                for (label jAngle = 0; jAngle < nAngle; jAngle++)
                {
                    label sweepRayID = jAngle + iBand * nAngle;
                    vector sweepDir = dom.IRay(sweepRayID).d();
                    vector sweepdAve = dom.IRay(sweepRayID).dAve();

                    scalar cosB = surfNorm& sweepDir;

                    if (cosB > 0.0)
                    {
                        vector reflecIncidentDir = sweepDir - 2 * cosB * surfNorm;
                        label reflecIncidentRay = -1;
                        dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                        const scalarField& reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];
                        diffusive = diffusive + reflecFace[faceI] * mag (surfNorm& sweepdAve);
                    }
                }
            }

            for (label i = 1; i <= npTheta; i++)
            {
                scalar pxRayTheta = bdRayTheta - 0.5 * deltaTheta + 0.5 * (2 * i - 1) * deltaTheta / npTheta;

                for (label j = 1; j <= npPhi; j++)
                {
                    scalar pxRayPhi = bdRayPhi - 0.5 * deltaPhi + 0.5 * (2 * j - 1) * deltaPhi / npPhi;
                    scalar sinTheta = Foam::sin(pxRayTheta);
                    scalar cosTheta = Foam::cos(pxRayTheta);
                    scalar sinPhi = Foam::sin(pxRayPhi);
                    scalar cosPhi = Foam::cos(pxRayPhi);
                    vector pixelDir = vector(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
                    scalar cosB = pixelDir& surfNorm;

                    if (cosB > 0.0)
                    {
                        scalar pixelOmega = 2.0 * sinTheta * Foam::sin(deltaTheta / 2.0/ npTheta) * deltaPhi / npPhi;
                        vector reflecIncidentDir = pixelDir - 2 * cosB * surfNorm;
                        label reflecIncidentRay = -1;
                        dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                        const scalarField& reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];
                        spectacular = spectacular + reflecFace[faceI] * pixelOmega;
                    }
                }
            }

            refValue()[faceI] = reflectionCoef_ * 
            (diffuseFraction_ * diffusive/ pi / 2 + (1.0 - diffuseFraction_) * spectacular / bdOmega);
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
        }
        else
        {
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0;
        }

    }
    
    // Restore tag
    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs();
}

void reflectiveSurfaceMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("reflectionCoef ") << reflectionCoef_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction ") << diffuseFraction_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    reflectiveSurfaceMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam


// ************************************************************************* //
