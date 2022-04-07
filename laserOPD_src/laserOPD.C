/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "laserOPD.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "interpolationCell.H"
#include "interpolationCellPoint.H"
#include "Random.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(laserOPD, 0);
        addToRadiationRunTimeSelectionTables(laserOPD);
    }

    defineTemplateTypeNameAndDebugWithName
    (
        Cloud<laserParticle>,
        "laserCloud",
        0
    );

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::laserOPD::initialise()
{
    // Initialise the laser particles
    laserCloud_.clear();

    if (mesh_.nGeometricD() == 3)
    {
        const vector lDir = normalised(laserDirection_);
        laserMode_.createCloud(laserCloud_, *this, maxTrackLength_, focalLaserPosition_, lDir);
    }
    else
    {
        FatalErrorInFunction
            << "Current functionality limited to 3-D cases"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::laserOPD::laserOPD(const volScalarField& rho)
:
    radiationModel(typeName, rho),
    laserCloud_(mesh_, "laserCloud", IDLList<laserParticle>()),
    // maxTrackLength_(get<scalar>("maxTrackLength")),
    // densityRef_(get<scalar>("densityRef")),
    // GladstoneDale_(get<scalar>("GladstoneDale")),
    // focalLaserPosition_(get<point>("focalLaserPosition")),
    // laserDirection_(get<vector>("laserDirection")),
    maxTrackLength_(readScalar(this->lookup("maxTrackLength"))),
    densityRef_(readScalar(this->lookup("densityRef"))),
    GladstoneDale_(readScalar(this->lookup("GladstoneDale"))),
    focalLaserPosition_(this->lookup("focalLaserPosition")),
    laserDirection_(this->lookup("laserDirection")),
    laserModePtr_(radiation::laserModel::New(*this, this->mesh_)),
    laserMode_(laserModePtr_()),
    GladstoneDaleField_
    (
        IOobject
        (
            "GD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimDensity, GladstoneDale_)
    ),
    n_
    (
        IOobject
        (
            "n",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    rho_(rho)
{
    initialise();
}


Foam::radiation::laserOPD::laserOPD
(
    const dictionary& dict,
    const volScalarField& rho
)
:
    radiationModel(typeName, dict, rho),
    laserCloud_(mesh_, "laserCloud", IDLList<laserParticle>()),
    // maxTrackLength_(get<scalar>("maxTrackLength")),
    // densityRef_(get<scalar>("densityRef")),
    // GladstoneDale_(get<scalar>("GladstoneDale")),
    // focalLaserPosition_(get<point>("focalLaserPosition")),
    // laserDirection_(get<vector>("laserDirection")),
    maxTrackLength_(readScalar(this->lookup("maxTrackLength"))),
    densityRef_(readScalar(this->lookup("densityRef"))),
    GladstoneDale_(readScalar(this->lookup("GladstoneDale"))),
    focalLaserPosition_(this->lookup("focalLaserPosition")),
    laserDirection_(this->lookup("laserDirection")),
    laserModePtr_(radiation::laserModel::New(*this, this->mesh_)),
    laserMode_(laserModePtr_()),
    GladstoneDaleField_
    (
        IOobject
        (
            "GD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimDensity, GladstoneDale_)
    ),
    n_
    (
        IOobject
        (
            "n",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    rho_(rho)
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::laserOPD::read()
{
    if (radiationModel::read())
    {
        return true;
    }

    return false;
}

Foam::label Foam::radiation::laserOPD::nBands() const
{
    return 1;
}


void Foam::radiation::laserOPD::calculate()
{

    // Clear and initialise the cloud
    // NOTE: Possible to reset original particles, but this requires
    // data transfer for the cloud in differet processors.
    initialise();

    // Reset the field
    n_ = GladstoneDaleField_*rho_ - GladstoneDale_*densityRef_;

    const interpolationCell<scalar> nInterp(n_);

    laserParticle::trackingData td
    (
        laserCloud_,
        nInterp
    );

    Info<< "Move particles..."
        << returnReduce(laserCloud_.size(), sumOp<label>()) << endl;

    laserCloud_.move(laserCloud_, td, mesh_.time().deltaTValue());

}


Foam::tmp<Foam::volScalarField> Foam::radiation::laserOPD::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        mesh_,
        dimensionedScalar(dimPower/dimVolume/pow4(dimTemperature), Zero)
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::laserOPD::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        mesh_,
        dimensionedScalar(dimPower/dimVolume, Zero)
    );
}


// ************************************************************************* //
