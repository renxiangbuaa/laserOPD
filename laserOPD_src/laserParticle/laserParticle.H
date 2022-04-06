/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd
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
    Foam::laserParticle

Description
    Discrete Transfer Radiation Model (laser) particle

SourceFiles
    laserParticle.H

\*---------------------------------------------------------------------------*/

#ifndef laserParticle_H
#define laserParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCell.H"
#include "volFieldsFwd.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class laserParticle;

// Forward declaration of friend functions
Ostream& operator<<(Ostream&, const laserParticle&);

/*---------------------------------------------------------------------------*\
                        Class laserParticle Declaration
\*---------------------------------------------------------------------------*/

class laserParticle
:
    public particle
{
    // Private data

        //- Relative coordinates (1 x y)
        vector xy_;

        //- Initial position
        point p0_;

        //- Target position
        point p1_;

        //- optical path difference [m]
        scalar OPL_;

        //- Area of radiation
        scalar dA_;

        //- 
        bool hitWall_;


public:

    friend class Cloud<laserParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCell<scalar>& nInterp_;


    public:

        // Constructors

            inline trackingData
            (
                Cloud<laserParticle>& spc,
                const interpolationCell<scalar>& nInterp
            );


        // Member functions

            inline const interpolationCell<scalar>& nInterp() const;
    };


    // Static Data Members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


        //- String representation of properties
        AddToPropertyList
        (
            particle,
            " xy"
          + " p0"
          + " p1"
          + " OPL"
          + " dA"
          + " hitWall_";
        );


    // Constructors

        //- Construct from components, with searching for tetFace and
        //  tetPt unless disabled by doCellFacePt = false.
        laserParticle
        (
            const polyMesh& mesh,
            const vector& xy,
            const vector& position,
            const vector& targetPosition,
            const label cellI,
            const scalar dA
        );

        //- Construct from components
        laserParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& xy,
            const vector& position,
            const vector& targetPosition,
            const scalar dA
        );

        //- Construct from Istream
        laserParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as copy
        laserParticle(const laserParticle& p);


        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<laserParticle> operator()(Istream& is) const
            {
                return autoPtr<laserParticle>
                (
                    new laserParticle(mesh_, is, true)
                );
            }
        };


    // Access

        //- Return const access to the relative coordinates (1 x y)
        inline const vector& xy() const;

        //- Return const access to the initial position
        inline const point& p0() const;

        //- Return const access to the target position
        inline const point& p1() const;

        //- Return const access to the current intensity
        inline scalar OPL() const;

        //- Return const access dA
        inline scalar dA() const;


    // Edit

        //- Return access to the target position
        inline point& p1();

        //- Return const access to the current intensity
        inline scalar& OPL();

        //- Return access to dA
        inline scalar& dA();


   // Tracking

        //- Move
        bool move(Cloud<laserParticle>& , trackingData&, const scalar);


    // Member Operators

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch
        (
            Cloud<laserParticle>&,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch
        (
            Cloud<laserParticle>&,
            trackingData& td
        );

        bool hitPatch
        (
            Cloud<laserParticle>&,
            trackingData& td
        );


            static void writeFields(const Cloud<laserParticle>& c);

    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const laserParticle& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "laserParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
