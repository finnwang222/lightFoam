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
#include "transExteriorSurfaceMixedFvPatchScalarField.H"
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

transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(0.0),
    nOwn_(0.0),
    I0_(0.0),
    nBands_(1),
    diffuseFraction_(0.0),
    beamWidthPhi_(0.0),
    beamWidthTheta_(0.0),
    beamDir_(Zero),
    beamNormToSurf_(false),
    initFlag_(0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const transExteriorSurfaceMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    diffuseFraction_(ptf.diffuseFraction_),
    beamWidthPhi_(ptf.beamWidthPhi_),
    beamWidthTheta_(ptf.beamWidthTheta_),
    beamDir_(ptf.beamDir_),
    beamNormToSurf_(ptf.beamNormToSurf_),
    initFlag_(ptf.initFlag_)
{}


transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(readScalar(dict.lookup("nNbg"))),     //scalarField("n1", dict),
    nOwn_(readScalar(dict.lookup("nOwn"))),     //scalarField("n2", dict)
    I0_(readScalar(dict.lookup("irradiation"))),
    nBands_(readLabel(dict.lookup("nBand"))),
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction"))),
    beamWidthPhi_(readScalar(dict.lookup("beamWidthPhi"))),
    beamWidthTheta_(readLabel(dict.lookup("beamWidthTheta"))),
    beamDir_(dict.lookup("beamDir")),
    initFlag_(0)
{
    lightBandDist_.setSize(nBands_);
    dict.lookup("lightBandDist") >> lightBandDist_;
    dict.lookup("beamNormToSurf") >> beamNormToSurf_;

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
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
        fvPatchScalarField::operator = (refValue());
    }
}

transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const transExteriorSurfaceMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    diffuseFraction_(ptf.diffuseFraction_),
    beamWidthPhi_(ptf.beamWidthPhi_),
    beamWidthTheta_(ptf.beamWidthTheta_),
    beamDir_(ptf.beamDir_),
    beamNormToSurf_(ptf.beamNormToSurf_),
    initFlag_(ptf.initFlag_)
{}

transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const transExteriorSurfaceMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    diffuseFraction_(ptf.diffuseFraction_),
    beamWidthPhi_(ptf.beamWidthPhi_),
    beamWidthTheta_(ptf.beamWidthTheta_),
    beamDir_(ptf.beamDir_),
    beamNormToSurf_(ptf.beamNormToSurf_),
    initFlag_(ptf.initFlag_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void transExteriorSurfaceMixedFvPatchScalarField::updateCoeffs()
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

    if(dom.nBand() == 0)
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

    label rayId = -1;
    dom.setRayId(internalField().name(), rayId);

    const label iBand = dom.IRay(rayId).iBand();
    const label nTheta = dom. nTheta();
    const label nPhi = dom.nPhi();
    const vector& bdRayDir = dom.IRay(rayId).d();
    const vector& bdRayDAve = dom.IRay(rayId).dAve();
    const scalar bdOmega = dom.IRay(rayId).omega();
    const scalar bdRayPhi = dom.IRay(rayId).phi();
    const scalar bdRayTheta = dom.IRay(rayId).theta();
    const scalar deltaPhi = pi/ (2.0 * nPhi);
    const scalar deltaTheta = pi / nTheta;

    const label nAngle = dom.nAngle();

    label npPhi = dom.NumPixelPhi();
    label npTheta = dom.NumPixelTheta();

    label numDivPhi = 1;
    label numDivTheta = 1;
    label rayJ = -1;
    scalar deltaDivTheta;
    scalar beamThetaStart;

    scalar nRatio = nOwn_ / nNbg_;

    dirToAngle(beamDir_, beamAnglePhi_, beamAngleTheta_);

    if(initFlag_ == 0)
    {
        initFlag_ = 1;
        spectacularRefraction.set(new scalarField(patch().size(), 0.0));
        diffusiveRefraction.set(new scalarField(patch().size(), 0.0));

        if (beamNormToSurf_)
        {
            forAll(Iw, faceI)
            {
                scalar T = 4 * nNbg_ * nOwn_ / (nNbg_ + nOwn_) / (nNbg_ + nOwn_);

                if(diffuseFraction_ > 0)
                {
                    diffusiveRefraction()[faceI] = diffusiveRefraction()[faceI] + I0_ * T * mag(n[faceI] & bdRayDAve);
                }

                dom.dirToRayId(-n[faceI], iBand, rayJ);
                if(rayJ == rayId)
                {
                    spectacularRefraction()[faceI] = spectacularRefraction()[faceI] + I0_ * T * bdOmega;
                }
            }
        }
        else
        {
            scalar bwP = beamWidthPhi_ / 180 * pi;
            scalar bwT = beamWidthTheta_ /180 * pi;

            if(internalField().mesh().nSolutionD() ==3)     //3D
            {
                if(bwP > deltaPhi)
                {
                    numDivPhi = npPhi * label (bwP / deltaPhi);
                }
                else
                {
                    if(bwP > deltaPhi/2)
                    {
                        numDivPhi = npPhi;
                    }
                    else
                    {
                        numDivPhi = 1;
                    }
                }

                if(bwT > deltaTheta)
                {
                    numDivTheta = npTheta * label (bwT / deltaTheta);
                }
                else
                {
                    if(bwT > deltaTheta/2)
                    {
                        numDivTheta = npTheta;
                    }
                    else
                    {
                        numDivTheta = 1;
                    }
                }
                deltaDivTheta = bwT / numDivTheta;
                beamThetaStart = beamAngleTheta_ - bwT / 2.0;
            }
            else
            {
                if(internalField().mesh().nSolutionD() == 2)    //2D
                {
                    if(bwP > deltaPhi)
                    {
                        numDivPhi = npPhi * label(bwP / deltaPhi);
                    }
                    else
                    {
                        if (bwP > deltaPhi /2)
                        {
                            numDivPhi = npPhi;
                        }
                        else
                        {
                            numDivPhi = 1;
                        }
                    }
                    numDivTheta = 1;
                    deltaDivTheta = pi;
                    beamThetaStart = 0;
                }
                else    //1D
                {
                    numDivPhi = 1;
                    numDivTheta = 1;
                    deltaDivTheta = pi;
                    beamThetaStart = 0;
                }
            }
            scalar deltaDivPhi = bwP / numDivPhi;
            scalar beamPhiStart = beamAnglePhi_ - bwP / 2.0;

            forAll(Iw, faceI)
            {
                for (label i =1; i <= numDivTheta; i++)
                {
                    scalar sweepDivTheta = beamPhiStart + (2 * i - 1) * deltaDivTheta / 2.0;
                    for (label j =1; j <= numDivPhi; j++)
                    {
                        scalar sweepDivPhi = beamThetaStart + (2 * j - 1) * deltaDivPhi / 2.0;
                        scalar sinTheta = Foam::sin(sweepDivTheta);
                        scalar cosTheta = Foam::cos(sweepDivTheta);
                        scalar sinPhi = Foam::sin(sweepDivPhi);
                        scalar cosPhi = Foam::cos(sweepDivPhi);
                        vector sweepDivDir = vector(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
                        vector surfNorm = -n[faceI];
                        scalar cosA = surfNorm& sweepDivDir;

                        if(cosA > 0.0 && cosA * cosA > 1-1/(nRatio * nRatio))   //direction out of wall
                        {
                            vector refracDir = (sweepDivDir - cosA * surfNorm) / nRatio
                            + Foam::sqrt(1 -1/(nRatio * nRatio) * (1 - cosA * cosA)) * surfNorm;

                            scalar cosB = surfNorm& refracDir;
                            scalar r1 = (nNbg_ * cosB - nOwn_ * cosA) / (nNbg_ * cosB + nOwn_ * cosA);
                            scalar r2 = (nNbg_ * cosA - nOwn_ * cosB) / (nNbg_ * cosA + nOwn_ * cosB);
                            scalar T = 1.0 - 0.5 * (r1 * r1 + r2 * r2);
                            scalar refracPhi;
                            scalar refracTheta;
                            dirToAngle(refracDir, refracPhi, refracTheta);
                            sinTheta = Foam::sin(sweepDivTheta);
                            cosTheta = Foam::cos(sweepDivTheta);
                            sinPhi = Foam::sin(sweepDivPhi);
                            cosPhi = Foam::cos(sweepDivPhi);

                            if(diffuseFraction_ > 0)
                            {
                                vector refracdAve = vector
                                (
                                    cosPhi * Foam::sin(0.5 * deltaDivPhi)
                                    * (deltaDivTheta - Foam::cos(2.0 * refracTheta))
                                    * Foam::sin(deltaTheta),
                                    sinPhi * Foam::sin(0.5 * deltaDivPhi)
                                    * (deltaDivTheta - Foam::cos(2.0 * refracTheta))
                                    * Foam::sin(deltaTheta),
                                    0.5 * deltaDivPhi
						            *Foam::sin(2.0*refracTheta)*Foam::sin(deltaDivTheta)
                                );
                                diffusiveRefraction()[faceI] = diffusiveRefraction()[faceI] + I0_ * T * mag(n[faceI]& refracdAve);
                            }

                            label refracRay = -1;
                            dom.dirToRayId(refracDir, iBand, refracRay);
                            if(refracRay == rayId)
                            {
                                scalar refracOmega = 2.0 * sinTheta * Foam::sin(deltaDivTheta / 2.0) * deltaDivPhi;
                                spectacularRefraction()[faceI] = spectacularRefraction()[faceI] + I0_ * T * refracOmega;
                            }
                        }
                    }
                }
            }
        }
    }

    scalar spectacularReflection = 0.0;
    scalar diffusiveReflection = 0.0;

    if (internalField().mesh().nSolutionD() == 2)
    {
        npTheta = 1;
    }
    if (internalField().mesh().nSolutionD() == 1)
    {
        npTheta = 1;
        npPhi =1;
    }

    forAll (Iw, faceI)
    {
        spectacularReflection = 0.0;
        diffusiveReflection = 0.0;
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
                    vector sweepdAve = dom.IRay(rayJ).dAve();

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
                            diffusiveReflection = diffusiveReflection +
                            reflecFace[faceI] * R * mag(surfNorm& sweepdAve);
                        }
                        else
                        {
                            diffusiveReflection = diffusiveReflection + 
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
                            spectacularReflection = spectacularReflection + 
                            reflecFace[faceI] * R * pixelOmega;
                        }
                        else
                        {
                            spectacularReflection = spectacularReflection + reflecFace[faceI] * pixelOmega;
                        }
                    }

                }
            }
            refValue()[faceI] = diffuseFraction_ * (diffusiveRefraction()[faceI] + diffusiveReflection) / pi / 2
            + (1.0 - diffuseFraction_) * (spectacularReflection + spectacularRefraction()[faceI])/ bdOmega;
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


void Foam::light::transExteriorSurfaceMixedFvPatchScalarField::dirToAngle
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
		if(dir.x() > 0 && dir.y() < 0 ) tPhi = tPhi + 2*pi;
	}
	else
	{
		if(dir.y() > 0 ) tPhi = pi/2.0;
		if(dir.y() < 0 ) tPhi = 3.0*pi/2.0;;
	}
	
}

void transExteriorSurfaceMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nNbg") << nNbg_ << token::END_STATEMENT << nl;
    os.writeKeyword("nOwn") << nOwn_ << token::END_STATEMENT <<nl;
    os.writeKeyword("I0") << I0_ << token::END_STATEMENT << nl;
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction") << diffuseFraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("beamWidthPhi_") << beamWidthPhi_ << token::END_STATEMENT << nl;
    os.writeKeyword("beamDir") << beamDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("beamNormToSurf") << beamNormToSurf_ << token::END_STATEMENT << nl;
                
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    transExteriorSurfaceMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam


// ************************************************************************* //
