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

\*---------------------------------------------------------------------------*/

#include "laserParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::laserParticle::propertyList_ =
    Foam::laserParticle::propertyList();

const std::size_t Foam::laserParticle::sizeofFields_
(
    sizeof(laserParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laserParticle::laserParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    xy_(Zero),
    p0_(Zero),
    p1_(Zero),
    OPL_(0),
    dA_(0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> xy_ >> p0_ >> p1_ >> OPL_ >> dA_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, xy_.data(), vector::nComponents);
            readRawScalar(is, p0_.data(), vector::nComponents);
            readRawScalar(is, p1_.data(), vector::nComponents);
            readRawScalar(is, &OPL_);
            readRawScalar(is, &dA_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&p0_), sizeofFields_);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::laserParticle::writeFields
(
    const Cloud<laserParticle>& c
)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<vector> xy(c.fieldIOobject("xy", IOobject::NO_READ), np);
    // IOField<point> p0(c.fieldIOobject("p0", IOobject::NO_READ), np);
    // IOField<point> p1(c.fieldIOobject("p1", IOobject::NO_READ), np);
    IOField<scalar> OPL(c.fieldIOobject("OPL", IOobject::NO_READ), np);
    IOField<scalar> dA(c.fieldIOobject("dA", IOobject::NO_READ), np);


    label i = 0;
    for (const laserParticle& p : c)
    {
        xy[i] = p.xy();
        // p0[i] = p.p0();
        // p1[i] = p.p1();
        OPL[i] = p.OPL();
        dA[i] = p.dA();
        ++i;
    }

    xy.write(np > 0);
    // p0.write(np > 0);
    // p1.write(np > 0);
    OPL.write(np > 0);
    dA.write(np > 0);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const laserParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.xy_
            << token::SPACE << p.p0_
            << token::SPACE << p.p1_
            << token::SPACE << p.OPL_
            << token::SPACE << p.dA_
            << token::SPACE << p.hitWall_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.p0_),
            laserParticle::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
