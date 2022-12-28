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
#include "transInteriorSurfaceCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "regionProperties.H"
#include "lightDOM.H"
#include "mappedPatchFieldBase.H"
#include "constants.H"
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

transInteriorSurfaceCoupledFvPatchScalarField::
transInteriorSurfaceCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(0.0),
    nOwn_(0.0),
    diffuseFraction_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


transInteriorSurfaceCoupledFvPatchScalarField::
transInteriorSurfaceCoupledFvPatchScalarField
(
    const transInteriorSurfaceCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    diffuseFraction_(ptf.diffuseFraction_)
{}


transInteriorSurfaceCoupledFvPatchScalarField::
transInteriorSurfaceCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(readScalar(dict.lookup("nNbg"))),     //scalarField("n1", dict),
    nOwn_(readScalar(dict.lookup("nOwn"))),     //scalarField("n2", dict)
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction")))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
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
    }
}


transInteriorSurfaceCoupledFvPatchScalarField::
transInteriorSurfaceCoupledFvPatchScalarField
(
    const transInteriorSurfaceCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    nNbg_(wtcsf.nNbg_),
    nOwn_(wtcsf.nOwn_),
    diffuseFraction_(wtcsf.diffuseFraction_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void transInteriorSurfaceCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
   const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch().patch());
   const label samplePatchI = mpp.samplePolyPatch().index();
   const polyMesh& nbrMesh = mpp.sampleMesh();
   const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];
   const mapDistribute& distMap = mpp.map();

    scalarField Ic(patchInternalField());
    scalarField& Iw = *this;

    const lightModel& light = db().lookupObject<lightModel>("lightProperties");
    const lightDOM& dom(refCast<const lightDOM>(light));

    label BDrayId = -1;
    dom.setRayId(internalField().name(), BDrayId);

    const label patchI = patch().index();
    vectorField n(patch().Sf() / patch().magSf());

    const label iBand = dom.IRay(BDrayId).iBand();
    const label nTheta = dom. nTheta();
    const label nPhi = dom.nPhi();
    const vector bdRayDir = dom.IRay(BDrayId).d();
    const scalar bdOmega = dom.IRay(BDrayId).omega();
    const scalar deltaPhi = pi/ (2.0 * nPhi);
    const scalar deltaTheta = pi / nTheta;

    const label nAngle = dom.nAngle();
    const scalar bdRayPhi = dom.IRay(BDrayId).phi();
    const scalar bdRayTheta = dom.IRay(BDrayId).theta();
    label npPhi = dom.NumPixelPhi();
    label npTheta = dom.NumPixelTheta();

    PtrList<scalarField> NbrRaySet(nAngle);

    for(label jAngle = 0; jAngle<nAngle; jAngle++)
    {
        NbrRaySet.set(jAngle, new scalarField(patch().size(), 0.0));
        label rayI = jAngle + iBand * nAngle;
        const transInteriorSurfaceCoupledFvPatchScalarField&
        nbrField = refCast<const transInteriorSurfaceCoupledFvPatchScalarField>
        (nbrPatch.lookupPatchField<volScalarField, scalar>(dom.IRay(rayI).I().name()));
        scalarField TcNbr(nbrField.patchInternalField());
        distMap.distribute(TcNbr);
        NbrRaySet[jAngle] = TcNbr;
    }

    scalar spectacular = 0.0;
    scalar diffusive = 0.0;

    if (internalField().mesh().nSolutionD() == 2)   //2D
    {
        npTheta = 1;
    }
    if (internalField().mesh().nSolutionD() == 1)    //1D
    {
        npTheta = 1;
        npPhi = 1;
    }

    scalar nRatio = nOwn_/ nNbg_;

    forAll (Iw, faceI)
    {
        spectacular = 0.0;
        diffusive = 0.0;
        vector surfNorm = -n[faceI];
        scalar cosA = surfNorm& bdRayDir;

        if (cosA > 0.0)     //dirction out of the wall
        {
            if (diffuseFraction_ > 0)
            {
                for (label jAngle = 0; jAngle < nAngle; jAngle++)
                {
                    label sweepRayId = jAngle + iBand * nAngle;
                    vector sweepDir = dom.IRay(sweepRayId).d();
                    vector sweepdAve = dom.IRay(sweepRayId).dAve();

                    scalar cosB = surfNorm& sweepDir;

                    if (cosB > 0.0)     //direction out of the wall
                    {
                        vector reflecIncidentDir = sweepDir =- 2 * cosB * surfNorm;
                        label reflecIncidentRay = -1;
                        dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                        const scalarField& reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];

                        if (cosB * cosB > 1 - 1/(nRatio * nRatio))
                        {
                            vector refracIncidentDir = (sweepDir - cosB * surfNorm) * nRatio
                            + Foam::sqrt(1 - (nRatio * nRatio) * (1 - cosB * cosB)) * surfNorm;
                            scalar cosA = surfNorm& reflecIncidentDir;
                            scalar r1 = (nNbg_ * cosB - nOwn_ * cosA) 
                                      / (nNbg_ * cosB + nOwn_ * cosA);
                            scalar r2 = (nNbg_ * cosA - nOwn_ * cosB)
                                      / (nNbg_ * cosA + nOwn_ * cosB);
                            scalar R = 0.5 * (r1 * r1 + r2 * r2);
                            label refracIncidentRay = -1;
                            dom.dirToRayId (refracIncidentDir, 0, refracIncidentRay);
                            diffusive = diffusive +
                            (NbrRaySet[refracIncidentRay][faceI] * (1-R) + reflecFace[faceI] * R) * mag(surfNorm & sweepdAve);
                        }
                        else
                        {
                            diffusive = diffusive + 
                            reflecFace[faceI] * mag(surfNorm & sweepdAve);
                        }   
                    }     

                }
            }

            for (label i = 1; i <= npTheta; i++)
            {
                scalar pxRayTheta = bdRayTheta - 0.5 * deltaTheta + 0.5 * (2 * i - 1) * deltaTheta / npTheta;
                for (label j = 1; j <= npPhi; j++)
                {
                    scalar pxRayPhi = bdRayPhi - 0.5 * deltaPhi + 0.5 *(2 * j - 1) * deltaPhi / npPhi;
                    scalar sinTheta = Foam::sin (pxRayTheta);
                    scalar cosTheta = Foam::cos (pxRayTheta);
                    scalar sinPhi = Foam::sin(pxRayPhi);
                    scalar cosPhi = Foam::cos (pxRayPhi);
                    vector pixelDir = vector (sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
                    scalar cosB = pixelDir& surfNorm;

                    if (cosB > 0.0)
                    {
                        scalar pixelOmega = 2.0 * sinTheta * Foam::sin(deltaTheta / 2.0/ npTheta) * deltaPhi/ npPhi;
                        vector reflecIncidentDir = pixelDir - 2 * cosB * surfNorm;
                        label reflecIncidentRay = -1;
                        dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                        const scalarField& reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];

                        if (cosB * cosB > 1 -1 / (nRatio * nRatio))
                        {
                            vector refracIncidentDir = (pixelDir - cosB * surfNorm) * nRatio
                            + Foam::sqrt(1 - (nRatio * nRatio) * (1 - cosB * cosB)) * surfNorm;
                            scalar cosA = surfNorm& pixelDir;
                            scalar r1 = (nNbg_ * cosB - nOwn_ * cosA)
                                      / (nNbg_ * cosB + nOwn_ * cosA);
                            scalar r2 = (nNbg_ * cosA - nOwn_ * cosB)
                                      / (nNbg_ * cosA + nOwn_ * cosB);
                            scalar R = 0.5 * (r1 * r1 + r2 * r2);
                            label refracIncidentRay = -1;
                            dom.dirToRayId (refracIncidentDir, 0, refracIncidentRay);
                            spectacular = spectacular + 
                            (NbrRaySet[refracIncidentRay][faceI] * (1-R) + reflecFace[faceI] * R) *pixelOmega;
                        }
                        else
                        {
                            spectacular = spectacular + reflecFace[faceI] * pixelOmega;
                        }
                    }

                }
            }
            refValue()[faceI] = diffuseFraction_ * diffusive / pi / 2
            + (1.0 - diffuseFraction_) * spectacular/ bdOmega;
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


void transInteriorSurfaceCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nNbg") << nNbg_ << token::END_STATEMENT << nl;
    os.writeKeyword("nOwn") << nOwn_ << token::END_STATEMENT <<nl;
    os.writeKeyword("diffuseFraction") << diffuseFraction_ << token::END_STATEMENT <<nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    transInteriorSurfaceCoupledFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam


// ************************************************************************* //
