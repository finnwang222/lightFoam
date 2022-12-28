/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "lightDOM.H"
#include "constants.H"
#include "unitConversion.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(lightDOM, 0);
        addToRunTimeSelectionTable
        (
            lightModel,
            lightDOM,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::lightDOM::lightDOM(const volScalarField& intensity)
:
    lightModel(typeName, intensity),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLuminousIntensity, 0.0)
    ),
    diffusionScatter_
    (
        IOobject
        (
            "diffusionScatter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLuminousIntensity, 0.0)
    ),
    nTheta_(coeffs_.get<label>("nTheta")),
    nPhi_(coeffs_.get<label>("nPhi")),
    nAngle_(0),
    nRay_(0),
    nBand_(coeffs_.getOrDefault<label>("nBand", 1)),
    NumPixelPhi_(coeffs_.getOrDefault<label>("NumpixelPhi", 1)),
    NumPixelTheta_(coeffs_.getOrDefault<label>("NumpixelTheta", 1)),
    GLambda_(nBand_),
    IRay_(0),
    tolerance_
    (
        coeffs_.getOrDefaultCompat<scalar>
        (
            "tolerance",
            {{"convergence", 1712}},
            0
        )
    ),
    maxIter_(coeffs_.getOrDefault<label>("maxIter", 50))
{
    Info << "lightDOM number of Bands" << nBand_ << endl;
    // 3D
    
    if (mesh_.nSolutionD() == 3)
    {
        nAngle_ = 4*nPhi_*nTheta_;
        nRay_ = nAngle_*nBand_;
        IRay_.setSize(nRay_);

        deltaPhi = pi/(2.0*nPhi_);
        deltaTheta = pi/nTheta_;

        label i = 0;
        for (label iBand = 0; iBand < nBand_; iBand++)
        {
            for (label n = 1; n <= nTheta_; n++)
            {
                for (label m = 1; m <= 4*nPhi_; m++)
                {
                    label iAngle = m-1 + (n-1)*4*nPhi_;
                    scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
                    scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                    IRay_.set
                    (
                        i,
                        new lightIntensityRay
                        (
                            *this,
                            mesh_,
                            iBand,
                            iAngle,
                            phii,
                            thetai,
                            deltaPhi,
                            deltaTheta
                        )
                    );
                    i++;
                }
            }
        }
        
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        const scalar thetai = piByTwo;
        deltaTheta = pi;
        nAngle_ = 4*nPhi_;
        nRay_ = nAngle_*nBand_;
        IRay_.setSize(nRay_);
        deltaPhi = pi/(2.0*nPhi_);
        label i = 0;
        for (label iBand = 0; iBand < nBand_; iBand++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                label iAngle = m-1;
                const scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new lightIntensityRay
                    (
                        *this,
                        mesh_,
                        iBand,
                        iAngle,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta
                    )
                );
                i++;
            }
        }
    }
    // 1D
    else
    {
        const scalar thetai = piByTwo;
        deltaTheta = pi;
        nAngle_ = 2;
        nRay_ = nAngle_*nBand_;
        IRay_.setSize(nRay_);
        deltaPhi = pi;
        label i = 0;
        for (label iBand = 0; iBand < nBand_; iBand++)
        {
             for (label m = 1; m <= 2; m++)
            {
                label iAngle = m-1;
                const scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new lightIntensityRay
                    (
                        *this,
                        mesh_,
                        iAngle,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        i
                    )
                );
            i++;
            }
        }
       
    }

    // Construct intensity field for each wavelength
    forAll(GLambda_, iBand)
    {
        GLambda_.set
        (
            iBand,
            new volScalarField
            (
                IOobject
                (
                    "GLambda_" + Foam::name(iBand) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                G_
            )
        );
    }

    sBand_.setSize(nBand_);
    aBand_.setSize(nBand_);
    // Construct absorption field for each wavelength
    forAll(aBand_, iBand)
    {
        aBand_[iBand] = extinction_->a(iBand);
        sBand_[iBand] = extinction_->s(iBand);
    }


    pf0_.setSize(nBand_);
    forAll(pf0_, iBand)
    {
        if (inScatter_->inScatter())
        {
            pf0_[iBand] = inScatter_->correct(1.0, iBand);
        }
        else
        {
            pf0_[iBand] = 0.0;
        }
    }

    Info<< "lightDOM : Allocated " << IRay_.size() << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::lightDOM::~lightDOM()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::light::lightDOM::read()
{
    if (lightModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresentCompat
        (
            "tolerance", {{"convergence", 1712}}, tolerance_
        );
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }

    return false;
}


void Foam::light::lightDOM::calculate()
{
    scalar maxResidual = 0.0;
    label rayJ;
    scalar rayCos;
    label iBand = 0;
    label iAngle = 0;
    label radIter = 0;
    scalar maxBandResidual = 0.0;

    do
    {
        

        radIter++;
        maxResidual = 0.0;
        forAll(IRay_, rayI)
        {
            iBand = IRay_[rayI].iBand();
            iAngle = IRay_[rayI].iAngle();

            Info<< "Light solver iter: " << radIter <<" iBand: "<<
            iBand << " iAngle: " << iAngle << endl;

            if (inScatter_->inScatter())
            {
                diffusionScatter_ = dimensionedScalar("diffusionScatter", dimLuminousIntensity, 0.0);

                for(label jAngle = 0; jAngle < nAngle_; jAngle++)
                {
                    rayJ = jAngle + iBand * nAngle_;
                    if(rayI != rayJ)
                    {
                        rayCos = IRay_[rayI].d() & IRay_[rayJ].d();

                        if(rayCos > 0)
                        {
                            diffusionScatter_ = diffusionScatter_ + IRay_[rayJ].I() * 
                            inScatter_->correct(rayCos, iBand) * IRay_[rayJ].omega()/2.0;
                        }
                    }
                }
            }
            IRay_[rayI].updateBoundary();
            maxBandResidual = IRay_[rayI].correct();
            maxResidual = max(maxBandResidual, maxResidual);

        }

    } while (maxResidual > tolerance_ && radIter < maxIter_);

    updateG();
}


void Foam::light::lightDOM::updateG()
{
    G_ = dimensionedScalar("zero", dimLuminousIntensity, 0.0);
    label rayI;

    forAll(GLambda_, iBand)
    {
        GLambda_[iBand] = dimensionedScalar("zero", dimLuminousIntensity, 0.0);

        for(label iAngle = 0; iAngle < nAngle_; iAngle++)
        {
            rayI = iAngle + iBand * nAngle_;
            GLambda_[iBand] += IRay_[rayI].I();
        }
        G_ += GLambda_[iBand];
    }
}


void Foam::light::lightDOM::setRayId
(
    const word& name,
    label& rayId
) const
{
    // Assuming name is in the form: CHARS_iBand_iAngle

    size_type i1 = name.find_first_of('_');
    size_type i2 = name.find_last_of('_');

    label ib = readLabel(IStringStream(name.substr(i1+1, i2-i1-1))());
    label ia = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());

    rayId = nAngle_ * ib +ia;
}

void Foam::light::lightDOM::dirToRayId
(
    const vector& dir,
    const label&  iBand,
    label& rayId
) const
{
    scalar tTheta = Foam::acos(dir.z() / mag(dir));
    scalar tPhi;
    if (dir.x() != 0)
    {
        tPhi = Foam::atan(dir.y() / dir.x());
        if (dir.x() < 0 && dir.y() >0)
        {
            tPhi = tPhi + pi;
        }
        if (dir.x() < 0 && dir.y() <0)
        {
            tPhi = tPhi + pi;
        }
         if (dir.x() > 0 && dir.y() <0)
        {
            tPhi = tPhi + 2 * pi;
        }
    }
    else
    {
        if(dir.y() > 0)
        {
            tPhi = pi / 2.0;
        }
        if(dir.y() < 0)
        {
            tPhi = 3 * pi / 2.0;
        }
    }
    label iPhi = label(tPhi / deltaPhi);
    label iTheta = label (tTheta / deltaTheta);
    rayId = nAngle_ * iBand + iTheta * 4 * nPhi_ + iPhi;
}

// ************************************************************************* //
