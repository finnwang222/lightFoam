/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd
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

#include "lightIntensityRay.H"
#include "fvm.H"
#include "lightDOM.H"
#include "constants.H"

using namespace Foam::constant;

Foam::label
Foam::light::lightIntensityRay::rayId(0);

const Foam::word
Foam::light::lightIntensityRay::intensityPrefix("I");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::lightIntensityRay::lightIntensityRay
(
    const lightDOM& dom,
    const fvMesh& mesh,
    const label iBand,
    const label iAngle,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta
)
:
    dom_(dom),
    mesh_(mesh),
    d_(Zero),
    dAve_(Zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    iBand_(iBand),
    iAngle_(iAngle)
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    d_ = vector(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
    dAve_ = vector
    (
        cosPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        0.5*deltaPhi*Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
    );

    autoPtr<volScalarField> IDefaultPtr;

    IOobject IHeader
    (
        intensityPrefix + "_" + name(iBand) + "_" + name(iAngle),
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    // Check if field exists and can be read
    if (IHeader.typeHeaderOk<volScalarField>(true))
    {
        I_.set
        (
            new volScalarField(IHeader, mesh_)
        );
    }
    else
    {
        // Demand driven load the IDefault field
        if (!IDefaultPtr)
        {
            IDefaultPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "IDefault",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_
                )
            );
        }

        // Reset the MUST_READ flag
        IOobject noReadHeader(IHeader);
        noReadHeader.readOpt(IOobject::NO_READ);

        I_.set
        (
            new volScalarField(noReadHeader, IDefaultPtr())
        );
    }
    rayId++;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::lightIntensityRay::~lightIntensityRay()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::light::lightIntensityRay::correct()
{
    scalar maxResidual = -GREAT;
    scalar eqnResidual;

    dimensionedScalar k("k", dimless/dimLength, dom_.a(iBand_) + dom_.s(iBand_) 
                        * (1.0 - dom_.pf0(iBand_) * omega_/2));
    dimensionedScalar s("s", dimless/dimLength, dom_.s(iBand_));

    const volScalarField& ds = dom_.diffusionScatter();
    
    const surfaceScalarField Ji(dAve_ & mesh_.Sf());

    fvScalarMatrix IiEq
    (
        fvm::div(Ji, I_(), "div(Ji,Ii_h)")
        + fvm::Sp(k*omega_, I_())
        == s * ds * omega_
    );

    IiEq.relax();

    eqnResidual = solve(IiEq, mesh_.solver("Ii")).initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    return maxResidual;
}


void Foam::light::lightIntensityRay::updateBoundary()
{
    I_->correctBoundaryConditions();
}


// ************************************************************************* //
