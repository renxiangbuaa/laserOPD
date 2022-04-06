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

#include "circularLaser.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(circularLaser, 0);
        addToRunTimeSelectionTable(laserModel, circularLaser, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::circularLaser::circularLaser
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    laserModel(dict, mesh),
    ndTheta_(dict.optionalSubDict(type() + "Coeffs").get<label>("nTheta")),
    ndr_(dict.optionalSubDict(type() + "Coeffs").get<label>("nr")),
    radius_(dict.optionalSubDict(type() + "Coeffs").get<scalar>("radius"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::circularLaser::createCloud(
    Cloud<laserParticle>& laserCloud,
    const dictionary& bd,
    const scalar maxTrackLength,
    const point lPosition,
    const vector lDir
)
{
    // Find a vector on the area plane. Normal to laser direction
    vector rArea = Zero;
    scalar magr = 0.0;

    {
        Random rnd(1234);

        while (magr < VSMALL)
        {
            vector v = rnd.sample01<vector>();
            rArea = v - (v & lDir)*lDir;
            magr = mag(rArea);
        }
    }
    rArea.normalise();

    scalar dr = radius_/ndr_;
    scalar dTheta = mathematical::twoPi/ndTheta_;
    

    // Count the number of missed positions
    label nMissed = 0;

    // Target position
    point p1 = vector::zero;

    // Seed laser particles
    point p0 = lPosition;

    for (label ri = 0; ri < ndr_; ri++)
    {
        scalar r1 = SMALL + dr*ri;

        scalar r2 = r1 + dr;

        scalar rP = ((r1 + r2)/2);

        // local radius on disk
        vector localR = rP*rArea;

        scalar theta0 = 0.0;//dTheta/2.0;
        for (label thetai = 0; thetai < ndTheta_; thetai++)
        {
            scalar theta1 = theta0 + SMALL  + dTheta*thetai;

            scalar theta2 = theta1 + dTheta;

            scalar thetaP = (theta1 + theta2)/2.0;

            quaternion Q(lDir, thetaP);

            // Initial position on disk
            vector initialPos = (Q.R() & localR);

            // Initial position
            p0 = lPosition + initialPos;

            // calculate target point using new deviation rl
            p1 = p0 + (maxTrackLength*lDir);

            scalar dAi = (sqr(r2) - sqr(r1))*(theta2 - theta1)/2.0;

            label cellI = this->mesh_.findCell(p0);

            vector xy(1,rP*cos(thetaP),rP*sin(thetaP));

            if (cellI != -1)
            {
                // Create a new particle
                laserParticle* pPtr =
                    new laserParticle(this->mesh_, xy, p0, p1, cellI, dAi);

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


// ************************************************************************* //
