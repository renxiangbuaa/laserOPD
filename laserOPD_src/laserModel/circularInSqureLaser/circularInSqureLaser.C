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

\*---------------------------------------------------------------------------*/

#include "circularInSqureLaser.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(circularInSqureLaser, 0);
        addToRunTimeSelectionTable(laserModel, circularInSqureLaser, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::circularInSqureLaser::circularInSqureLaser
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    laserModel(dict, mesh),
    nLength_(readLabel(dict.optionalSubDict(type() + "Coeffs").lookup("nLength"))),
    length_(readScalar(dict.optionalSubDict(type() + "Coeffs").lookup("length"))),
    directionX_(dict.optionalSubDict(type() + "Coeffs").lookup("directionX")),
    radius_(readScalar(dict.optionalSubDict(type() + "Coeffs").lookup("radius")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::circularInSqureLaser::createCloud(
    Cloud<laserParticle>& laserCloud,
    const dictionary& bd,
    const scalar maxTrackLength,
    const point lPosition,
    const vector lDir
)
{
    scalar dx = length_/nLength_;
    scalar dAi = sqr(dx);
    const vector xDir = normalised(directionX_);
    const vector yDir = xDir ^ lDir;
    // if((xDir & lDir) == 0)
    if (fabs((xDir^yDir)&lDir) == 1)
    {
        // Count the number of missed positions
        label nMissed = 0;

        // Target position
        point p1 = vector::zero;
        point p0 = vector::zero;

        for (label xi = 0; xi < nLength_; xi++)
        {
            scalar xP = dx*(xi+0.5) - 0.5*length_;

            vector localX = xP*xDir;

            for (label yi = 0; yi < nLength_; yi++)
            {

                scalar yP = dx*(yi+0.5) - 0.5*length_ ;

                vector localY = yP*yDir;

                p0 = lPosition + (localX + localY);

                p1 = p0 + (maxTrackLength*lDir);

                label cellI = mesh_.findCell(p0);

                vector xy(1, xP, yP);

                if (cellI != -1)
                {
                    // Create a new particle
                    scalar aa = dAi;
                    if((sqr(xP)+sqr(yP)) > sqr(radius_))
                    {
                        aa=0;
                    }

                    laserParticle* pPtr =
                        new laserParticle(mesh_, xy, p0, p1, cellI, aa);

                    // Add to cloud
                    laserCloud.addParticle(pPtr);
                }

                if (returnReduce(cellI, maxOp<label>()) == -1)
                {
                    if (++nMissed <= 10)
                    {
                        WarningInFunction
                            << "Cannot find owner cell for focalPoint at "
                            << p0 << endl;
                    }
                }
            }
        }

        if (nMissed)
        {
            Info<< "Seeding missed " << nMissed << " locations" << endl;
        }
        
    }
    else
    {
        Info << "directionX and laserDirection are not perpendicular." << endl;
    }



}


// ************************************************************************* //
