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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::laserParticle::trackingData::trackingData
(
    Cloud<laserParticle>& spc,
    const interpolationCell<scalar>& nInterp
)
:
    particle::trackingData(spc),
    nInterp_(nInterp)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCell<Foam::scalar>&
Foam::laserParticle::trackingData::nInterp() const
{
    return nInterp_;
}


inline const Foam::vector& Foam::laserParticle::xy() const
{
    return xy_;
}


inline const Foam::point& Foam::laserParticle::p0() const
{
    return p0_;
}


inline const Foam::point& Foam::laserParticle::p1() const
{
    return p1_;
}


inline Foam::scalar Foam::laserParticle::OPL() const
{
    return OPL_;
}


inline Foam::scalar Foam::laserParticle::dA() const
{
    return dA_;
}



inline Foam::scalar& Foam::laserParticle::dA()
{
    return dA_;
}


inline Foam::point& Foam::laserParticle::p1()
{
    return p1_;
}


inline Foam::scalar& Foam::laserParticle::OPL()
{
    return OPL_;
}

// ************************************************************************* //
