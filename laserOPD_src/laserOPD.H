/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::radiation::laserOPD

Group
    grpRadiationModels

Description
    Discrete Tray Radiation Method for collimated radiation flux.
    At the moment the particles are injected on a 2D disk in polar
    coordinates providing starting and finish locations.

    The ray is considered planar (no area divergence is considered). This
    model was developed to represent a collimated laser beam.


SourceFiles
    laserOPD.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModellaserOPD_H
#define radiationModellaserOPD_H

#include "radiationModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Cloud.H"
#include "laserParticle.H"
#include "laserModel.H"
#include "Function1.H"
#include "interpolation2DTable.H"
#include "labelField.H"
// #include "phasePairKey.H"
// #include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class laserOPD Declaration
\*---------------------------------------------------------------------------*/

class laserOPD
:
    public radiationModel
{

private:

    // Private data

        //- laser particle cloud
        Cloud<laserParticle> laserCloud_;

        //- Maximum tracking length for particles
        scalar maxTrackLength_;

        //- Reference density
        scalar densityRef_;

        //- Gladstone-Dale constant
        scalar GladstoneDale_;

        //- Focal laser position
        point focalLaserPosition_;

        //- Direction
        vector laserDirection_;


        // laser beam mode
        autoPtr<laserModel> laserModePtr_;
        laserModel& laserMode_;


        // Fields

            //- Gladstone-Dale constant
            volScalarField GladstoneDaleField_;

            //- refractive index
            volScalarField n_;

            //- 密度基准定度法
            const volScalarField& rho_;


    // Private Member Functions

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("laserOPD");


    // Constructors

        //- Construct from components
        laserOPD(const volScalarField& rho);

        //- Construct from components
        laserOPD(const dictionary& dict, const volScalarField& rho);


    //- Destructor
    virtual ~laserOPD() = default;


    // Member functions

    // Edit

        //- Solve radiation equation(s)
        void calculate();

        //- Read radiation properties dictionary
        bool read();

        //- Number of bands for this radiation model
        virtual label nBands() const;


    // Access

        //- Source term component (for power of T^4)
        virtual tmp<volScalarField> Rp() const;

        //- Source term component (constant)
        virtual tmp<DimensionedField<scalar, volMesh>> Ru() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
