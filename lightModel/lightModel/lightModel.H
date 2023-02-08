/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Namespace
    Foam::light

Description
    Namespace for light modelling

Class
    Foam::light::lightModel

Description
    Top level model for light modelling

SourceFiles
    lightModel.C
    lightModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef lightModel_H
#define lightModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatricesFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

// Forward declaration of classes
class extinctionModel;
class inScatterModel;

/*---------------------------------------------------------------------------*\
                       Class lightModel Declaration
\*---------------------------------------------------------------------------*/

class lightModel
:
    public IOdictionary
{
protected:



    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- light model on/off flag
        Switch light_;

        //- light model dictionary
        dictionary coeffs_;

        //- light solver frequency - number flow solver iterations per
        //  light solver iteration
        label solverFreq_;


        // References to the light sub-models

            //- extinction model
            autoPtr<extinctionModel> extinction_;

            //- inScatter model
            autoPtr<inScatterModel> inScatter_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- No copy construct
        lightModel(const lightModel&) = delete;

        //- No copy assignment
        void operator=(const lightModel&) = delete;


public:

    //- Runtime type information
    TypeName("lightModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            lightModel,
            dictionary,
            (
                const volScalarField& intensity
            ),
            (intensity)
        );


    // Constructors

        //- Null constructor
        lightModel(const volScalarField& intensity);

        //- Construct from components
        lightModel(const word& type, const volScalarField& intensity);



    // Selectors

        //- Return a reference to the selected light model
        static autoPtr<lightModel> New(const volScalarField& intensity);


    //- Destructor
    virtual ~lightModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve light equation(s)
            virtual void calculate() = 0;

            //- Read lightProperties dictionary
            virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToLightRunTimeSelectionTables(model)                            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        lightModel,                                                        \
        model,                                                                 \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        lightModel,                                                        \
        model,                                                                 \
        intensity                                                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
