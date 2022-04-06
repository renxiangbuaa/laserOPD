/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::radiation::circularLaser

Description
    Dummy scatter model for 'none'

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_circularLaser_H
#define radiation_circularLaser_H

#include "laserModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                         Class circularLaser Declaration
\*---------------------------------------------------------------------------*/

class circularLaser
:
    public laserModel
{
    // Private data

        //- Number of theta angles
        label ndTheta_;

        //- Number of radials
        label ndr_;

        //- Focal lase radius
        scalar radius_;

public:

    //- Runtime type information
    TypeName("circularLaser");

    //- Construct from components
    circularLaser(const dictionary& dict, const fvMesh& mesh);

    //- Destructor
    virtual ~circularLaser() = default;


    // Member Functions

        //- Read the laserMode dictionary
        virtual void createCloud(
            Cloud<laserParticle>&,
            const dictionary&,
            const scalar,
            const point,
            const vector
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
