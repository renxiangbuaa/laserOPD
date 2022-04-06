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

\*---------------------------------------------------------------------------*/

#include "laserParticle.H"
#include "constants.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laserParticle::laserParticle
(
    const polyMesh& mesh,
    const vector& xy,
    const vector& position,
    const vector& targetPosition,
    const label cellI,
    const scalar dA
)
:
    particle(mesh, position, cellI),
    xy_(xy),
    p0_(position),
    p1_(targetPosition),
    OPL_(0.),
    dA_(dA),
    hitWall_(false)
{}


Foam::laserParticle::laserParticle
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
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    xy_(xy),
    p0_(position),
    p1_(targetPosition),
    OPL_(0.),
    dA_(dA),
    hitWall_(false)
{}


Foam::laserParticle::laserParticle(const laserParticle& p)
:
    particle(p),
    xy_(p.xy_),
    p0_(p.p0_),
    p1_(p.p1_),
    OPL_(p.OPL_),
    dA_(p.dA_),
    hitWall_(p.hitWall_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::laserParticle::move
(
    Cloud<laserParticle>& spc,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    do
    {
        //Cache old data of particle to use for reflected particle
        const point pos0 = position();

        // Info<< pos0 << endl;

        const label cell1 = cell();
        const tetIndices tetIs = this->currentTetIndices();

        scalar f = 1 - stepFraction();
        const vector s = p1() - p0() - deviationFromMeshCentre();
        trackToAndHitFace(f*s, f, spc, td);

        const point p1 = position();
        vector dsv = p1 - pos0;
        scalar ds = mag(dsv);

        //const label cell1 = cell();

        //NOTE:
        // Under the new barocentric tracking alghorithm the newly
        // inserted particles are tracked to the nearest cell centre first,
        // then, given the direction, to a face. In both occasions the first call
        // to trackToAndHitFace returns ds = 0. In this case we do an extra
        // call to trackToAndHitFace to start the tracking.
        // This is a temporary fix until the tracking can handle it.
        if (ds == 0)
        {
            trackToAndHitFace(f*s, f, spc, td);
            dsv = p1 - position();
            ds = mag(dsv);
        }

        scalar n = td.nInterp().interpolate(pos0, cell1);

        OPL_ += n*ds;

    }while (!hitWall_ && !td.switchProcessor && stepFraction() < 1);

    // Info<< OPL_ << endl; 

    return td.keepParticle;
}



void Foam::laserParticle::hitProcessorPatch
(
    Cloud<laserParticle>&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::laserParticle::hitWallPatch
(
    Cloud<laserParticle>&,
    trackingData& td
)
{
    hitWall_= true;
    // td.keepParticle = false;
}


bool Foam::laserParticle::hitPatch
(
    Cloud<laserParticle>&,
    trackingData& td
)
{
    return false;
}


// ************************************************************************* //
