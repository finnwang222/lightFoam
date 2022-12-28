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
#include "lampMixedFvPatchScalarField.H"
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

lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    I0_(0.0),
    nBands_(1),
    lampRadius_(0.0),
    reflectionOnSurface_(false),
    reflectionCoef_(0.0),
    diffuseFraction_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const lampMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    lampRadius_(ptf.lampRadius_),
    reflectionOnSurface_(ptf.reflectionOnSurface_),
    reflectionCoef_(ptf.reflectionCoef_),
    diffuseFraction_(ptf.diffuseFraction_),
    lightBandDist_(ptf.lightBandDist_)
{}


lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    I0_(readScalar(dict.lookup("irradiation"))),
    nBands_(readLabel(dict.lookup("nBands"))),
    lampRadius_(readScalar(dict.lookup("lampRadius"))),
    reflectionCoef_(readScalar(dict.lookup("reflectionCoef"))),
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction")))
{
    dict.lookup("reflectionOnSurface") >> reflectionOnSurface_;
    lightBandDist_.setSize(nBands_);
    dict.lookup("lightBandDist") >> lightBandDist_;

    if (dict.found("value"))
    {
        // Full restart
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;
        fvPatchScalarField::operator = (refValue());
    }
}

lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const lampMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    lampRadius_(ptf.lampRadius_),
    reflectionOnSurface_(ptf.reflectionOnSurface_),
    reflectionCoef_(ptf.reflectionCoef_),
    diffuseFraction_(ptf.diffuseFraction_),
    lightBandDist_(ptf.lightBandDist_)
{}

lampMixedFvPatchScalarField::
lampMixedFvPatchScalarField
(
    const lampMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    lampRadius_(ptf.lampRadius_),
    reflectionOnSurface_(ptf.reflectionOnSurface_),
    reflectionCoef_(ptf.reflectionCoef_),
    diffuseFraction_(ptf.diffuseFraction_),
    lightBandDist_(ptf.lightBandDist_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void lampMixedFvPatchScalarField::updateCoeffs()
{
    if (this -> updated())
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

    if(dom.nBand() == 0 || internalField().mesh().nSolutionD() != 2)
    {
        FatalErrorIn
        (
            "Foam::light::"
            "wideBandDiffusiveRadiationMixedFvPatchScalarField::updateCoeffs"
        )
        <<"a non-grey bpundary condition is used with a grep"
        <<"absorption model" << nl << exit(FatalError);
    }

    const label patchI = patch().index();
    const vectorField n(patch().Sf() / patch().magSf());
    const scalarField sfSize = patch().magSf();

    label rayId = -1;
    dom.setRayId(internalField().name(), rayId);

    const label nAngle = dom.nAngle();
    const label iBand = dom.IRay(rayId).iBand();
    const vector& bdRayDir = dom.IRay(rayId).d();
    const scalar bdOmega = dom.IRay(rayId).omega();
    const scalar bdRayPhi = dom.IRay(rayId).phi();
    const scalar bdRayTheta = dom.IRay(rayId).theta();
    const label nTheta = dom. nTheta();
    
    const label nPhi = dom.nPhi();
    
    const scalar deltaPhi = pi/ (2.0 * nPhi);
    const scalar deltaTheta = pi / nTheta;

    label npPhi = dom.NumPixelPhi();
    label npTheta = dom.NumPixelTheta();

    const scalar nOwn_ = 1.00;
    const scalar nNbg_ = 1.5168;    //outside of air is glass
    const scalar nRatio = nOwn_ / nNbg_;
    const scalar height = 0.001;

    scalar spectacularReflection = 0.0;
    scalar reflectionRefraction = 0.0;
    scalar diffuseRefractionInside = 0.0;
    autoPtr<scalarField> diffusiveReflection;
    diffusiveReflection.set(new scalarField(patch().size(), 0.0));

    if(reflectionOnSurface_)
    {
        forAll (Iw, faceI)
        {
        
            vector surfNorm = -n[faceI];
        
            for (label jAngle = 0; jAngle < nAngle; jAngle++)
            {
                label sweepRayId = jAngle + iBand * nAngle;
                vector sweepDir = dom.IRay(sweepRayId).d();
                vector sweepdAve = dom.IRay(sweepRayId).dAve();
                scalar sweepdOmega = dom.IRay(sweepRayId).omega();

                scalar cosA = surfNorm& sweepDir;

                if (cosA < 0.0 && cosA >= -1)     //direction into of the wall
                {
                    const scalarField& reflecFace = dom.IRay(sweepRayId).I().boundaryField()[patchI];

                    scalar sinA = Foam::sqrt(1 - cosA * cosA);
                    scalar sinB = sinA * nRatio;
                    scalar cosB = -Foam::sqrt(1- sinB * sinB);
                        
                    scalar r1 = (nNbg_ * cosB - nOwn_ * cosA) 
                                / (nNbg_ * cosB + nOwn_ * cosA);
                    scalar r2 = (nNbg_ * cosA - nOwn_ * cosB)
                                / (nNbg_ * cosA + nOwn_ * cosB);
                    scalar R = 0.5 * (r1 * r1 + r2 * r2);
                            
                    diffusiveReflection()[faceI] = diffusiveReflection()[faceI] 
                                + reflecFace[faceI] * R * mag(n[faceI]& sweepdAve);
                    diffuseRefractionInside = diffuseRefractionInside 
                                + reflecFace[faceI] * (1.0 - R) * sfSize[faceI] *sweepdOmega;
                }        
            }       
        }
    }

        forAll(Iw, faceI)
        {
            vector surfNorm = -n[faceI];
            scalar cosA = surfNorm& bdRayDir;

            if(cosA > 0)
            {
                reflectionRefraction = 0.0;

                if(reflectionOnSurface_)
                {
                    spectacularReflection = 0.0;
                    scalar pxRayTheta = bdRayTheta;
                    for(label j =1; j <= npPhi; j++)
                    {
                        scalar pxRayPhi = bdRayPhi - 0.5 * deltaPhi + 0.5 * (2 * j - 1) * deltaPhi / npPhi;
                        scalar sinTheta = Foam::sin(pxRayTheta);
                        scalar cosTheta = Foam::cos(pxRayTheta);
                        scalar sinPhi = Foam::sin(pxRayPhi);
                        scalar cosPhi = Foam::cos(pxRayPhi);

                        vector pixelDir = vector(sinTheta * cosPhi, sinTheta *sinPhi, cosTheta);

                        scalar cosB = pixelDir& surfNorm;

                        if (cosB > 0.0)
                        {
                            scalar pixelOmega = 2.0 * sinTheta * Foam::sin(deltaTheta / 2.0/ npTheta) * deltaPhi/ npPhi;
                            vector reflecIncidentDir = pixelDir - 2 * cosB * surfNorm;
                            label reflecIncidentRay = -1;
                            dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                            const scalarField& reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];

                            if(cosB * cosB > 1-1/ (nRatio * nRatio))    //weather total reflection or not
                            {
                                vector refracIncidentDir = (pixelDir - cosB * surfNorm) * nRatio
                                    + Foam::sqrt(1-(nRatio * nRatio) * (1 - cosB * cosB)) * surfNorm;

                                scalar cosA = surfNorm& refracIncidentDir;     //mag(d);
                                scalar r1 = (nNbg_ * cosB - nOwn_ * cosA)
                                        / (nNbg_ * cosB + nOwn_ * cosA);
                                scalar r2 = (nNbg_ * cosA - nOwn_ * cosB)
                                        / (nNbg_ * cosA + nOwn_ * cosB);
                                scalar R = 0.5 * (r1 * r1 + r2 * r2);
                                spectacularReflection = spectacularReflection + reflecFace[faceI] * R * pixelOmega;
                            }
                            else
                            {
                                spectacularReflection = spectacularReflection +reflecFace[faceI] * pixelOmega;
                            }
                        }
                        reflectionRefraction = diffuseRefractionInside / (2 * pi * lampRadius_ * height) / 2 / pi * bdOmega
                                        + diffuseFraction_ * diffusiveReflection()[faceI] / pi /2
                                        + (1.0 - diffuseFraction_) * spectacularReflection / bdOmega;
                    }

                    refValue()[faceI] = reflectionCoef_ * reflectionRefraction + I0_ * lightBandDist_[iBand] / 2/ pi * bdOmega;
                    refGrad()[faceI] = 0.0;
                    valueFraction()[faceI] = 1.0;
                }
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


void Foam::light::lampMixedFvPatchScalarField::dirToAngle
(
    const vector& dir,
    scalar&	tPhi,
    scalar& tTheta
) const
{

	tTheta = Foam::acos(dir.z()/mag(dir));


	if(dir.x() != 0 )
	{ 
		tPhi = Foam::atan(dir.y()/dir.x());
		if(dir.x() < 0 && dir.y() >= 0 ) tPhi = tPhi + pi;
		if(dir.x() < 0 && dir.y() < 0 ) tPhi = tPhi + pi;
		if(dir.x() > 0 && dir.y() < 0 ) tPhi = tPhi + 2 * pi;
	}
	else
	{
		if(dir.y() > 0 ) tPhi = pi/ 2.0;
		if(dir.y() < 0 ) tPhi = 3.0 * pi/ 2.0;;
	}
	
}

void lampMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("irradiation") << I0_ << token::END_STATEMENT << nl;
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction") << diffuseFraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("lampRadius") << lampRadius_ << token::END_STATEMENT << nl;
    os.writeKeyword("reflectionOnSurface") << reflectionOnSurface_ << token::END_STATEMENT << nl;
    os.writeKeyword("reflectionCoef") << reflectionCoef_ << token::END_STATEMENT << nl;
    os.writeKeyword("lightBandDist") << lightBandDist_ << token::END_STATEMENT<< nl;
                
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    lampMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam


// ************************************************************************* //
