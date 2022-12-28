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

Class
    Foam::light::lightDOM

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media, not including scatter and reflective
    walls.

    Available absorption models:

    - constantAbsorptionEmission
    - greyMeanAbsoprtionEmission
    - wideBandAbsorptionEmission
    - multiBandAbsorptionEmission

    This model can handle non-grey participating media using
    multiBandAbsorptionEmission model. Accordingly the BC for rays should
    be wideBandDiffussive type

Usage
    \verbatim
        lightDOMCoeffs
        {
            nPhi        1;          // azimuthal angles in PI/2 on X-Y.
                                    //(from Y to X)
            nTheta      2;          // polar angles in PI (from Z to X-Y plane)
            tolerance   1e-3;       // convergence tolerance for light
                                    // iteration
        }

        solverFreq   1; // Number of flow iterations per radiation iteration
    \endverbatim

    The total number of solid angles is  4*nPhi*nTheta in 3-D.

    Operating modes:
    - 1-D:
      - ray directions are on X, Y or Z
      - \c nPhi and \c nTheta entries are ignored
      - \c meshOrientation vector can be used for any other 1-D direction.
    - 2-D:
      - ray directions are on X-Y, X-Z or Y-Z planes
      - only the \c nPhi entry is considered
      - \c meshOrientation vector can be used for non-aligned planes
        specifying the plane normal vector.
    - 3-D:
      - rays geberated in 3-D using the \c nPhi and \c nTheta entries
      - \c meshOrientation vector is not applicable.

    useSolarLoad calculates the primary and diffusive Sun fluxes on walls in
    addition to the RTE equations

    useExternalBeam add an external collimated beam to the domain. This option
    is not available if useSolarLoad is true.

    spectralDistribution is the energy spectral distribution of the collimated
    external beam.

SourceFiles
    lightDOM.C

\*---------------------------------------------------------------------------*/

#ifndef lightModellightDOM_H
#define lightModellightDOM_H

#include "lightIntensityRay.H"
#include "lightModel.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

/*---------------------------------------------------------------------------*\
                           Class lightDOM Declaration
\*---------------------------------------------------------------------------*/

class lightDOM
:
    public lightModel
{
    // Private Data

        //- Incident light  [W/m2]
        volScalarField G_;

        //boundary diffusion and inScattering light
        volScalarField diffusionScatter_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_ ;

        //- Number of solid angles
        label nAngle_;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nBand_;

        label NumPixelPhi_;

        label NumPixelTheta_;

        //- Wavelength total extinction coefficient [1/m]
        PtrList<volScalarField> GLambda_;

        DynamicList<scalar> aBand_;

        DynamicList<scalar> sBand_;

        DynamicList<scalar> pf0_;

        //- List of pointers to light intensity rays
        PtrList<lightIntensityRay> IRay_;

        //- Convergence tolerance
        scalar tolerance_;

        //- Maximum number of iterations
        label maxIter_;

        scalar deltaPhi;

        scalar deltaTheta;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- No copy construct
        lightDOM(const lightDOM&) = delete;

        //- No copy assignment
        void operator=(const lightDOM&) = delete;



public:

    //- Runtime type information
    TypeName("lightDOM");


    // Constructors

        //- Construct from components
        lightDOM(const volScalarField& intensity);

    //- Destructor
    virtual ~lightDOM();


    // Member Functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Update G and calculate total heat flux on boundary
            void updateG();

            //- Set the rayId and lambdaId from by decomposing an intensity
            //  field name
            void setRayId
            (
                const word& name,
                label& rayId
            ) const;

            void dirToRayId
            (
                const vector& dir,
                const label&  iBand,
                label& rayId
            ) const;


            //- Rotate rays according to Sun direction
            void updateRaysDir();

            //- Rotate rays spheric equator to sunDir
            void rotateInitialRays(const vector& sunDir);

            //- Align closest ray to sunDir
            void alignClosestRayToSun(const vector& sunDir);

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;


        // Access
            //- Ray intensity for rayI
            inline const lightIntensityRay& IRay(const label rayI) const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of angles
            inline label nAngle() const;

            //- Number of rays
            inline label nRay() const;

            //- Number of wavelength
            inline label nBand() const;

            inline label NumPixelPhi() const;

            inline label NumPixelTheta() const;

            inline scalar a(const label lambdaI) const;

            inline scalar s(const label lambdaI) const;

            inline scalar pf0(const label lambdaI) const;

            //- Const access to incident radiation field
            inline const volScalarField& G() const;

            inline const volScalarField& diffusionScatter() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lightDOMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
