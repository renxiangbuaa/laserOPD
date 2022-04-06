/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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

#include "laserOPDInfo.H"
#include "laserParticle.H"
#include "dictionary.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "tensor.H"
#include "QRMatrix.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(laserOPDInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        laserOPDInfo,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::laserOPDInfo::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Cloud information");
    writeCommented(os, "Time");
    writeTabbed(os, "meanOPD");
    writeTabbed(os, "coeffForTiltX");
    writeTabbed(os, "coeffForTiltY");
    writeTabbed(os, "rmsOPD");
    writeTabbed(os, "rmsOPDtilt");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::laserOPDInfo::laserOPDInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::laserOPDInfo::read(const dictionary& dict)
{
    if (regionFunctionObject::read(dict) && logFiles::read(dict))
    {
        logFiles::resetNames(dict.get<wordList>("clouds"));

        Info<< type() << " " << name() << ": ";
        if (writeToFile() && names().size())
        {
            Info<< "applying to clouds:" << nl;
            forAll(names(), cloudi)
            {
                Info<< "    " << names()[cloudi] << nl;
                writeFileHeader(files(cloudi));
            }
            Info<< endl;
        }
        else
        {
            Info<< "no clouds to be processed" << nl << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::laserOPDInfo::execute()
{
    return true;
}


bool Foam::functionObjects::laserOPDInfo::write()
{
    forAll(names(), cloudi)
    {
        const word& cloudName = names()[cloudi];

        const Cloud<laserParticle>& cloud =
            obr_.lookupObject<Cloud<laserParticle>>(cloudName);

        //- tip/tilt Coefficients
        tensor T(0, 0, 0, 0, 0, 0, 0, 0, 0);
        vector B(0, 0, 0);
        for (const laserParticle& p : cloud)
        {
            vector pp = p.xy();
            T += p.dA() * (pp * pp);
            B += (p.OPL()*p.dA()) * pp;
        }
        reduce(T,sumOp<tensor>());
        reduce(B,sumOp<vector>());
        scalar area=T.xx();  //T.xx() is The beam area
        T = T/area; 
        B = B/area;
        vector coeff=(inv(T) & B);

        Info << "The average optical path is: " << B.x() << endl;
        Info << "T: " << T << endl;
        Info << "B: " << B << endl;




        //- root mean square of OPD
        scalar rmsOPD=0;
        scalar rmsOPDtilt=0;
        for (const laserParticle& p : cloud)
        {
            rmsOPD += sqr(p.OPL()-coeff.x())*p.dA();
            rmsOPDtilt += sqr(p.OPL()-(coeff & p.xy()))*p.dA();
        }
        reduce(rmsOPD,sumOp<scalar>());
        reduce(rmsOPDtilt,sumOp<scalar>());
        rmsOPD = sqrt(rmsOPD/area);
        rmsOPDtilt = sqrt(rmsOPDtilt/area);
        
        Info << "The rmsOPD is: " << rmsOPD << endl;
        Info << "The rmsOPDtilt is: " << rmsOPDtilt << endl;

        // Log << type() << " " << name() <<  " write:" << nl
        //     << "    meanOPD    : " << meanOPD << nl
        //     << "    rmsOPD     : " << rmsOPD << nl
        //     << endl;




        if (writeToFile())
        {
            auto& os = files(cloudi);

            writeCurrentTime(os);
            os
                << token::TAB << coeff.x()
                << token::TAB << coeff.y()
                << token::TAB << coeff.z()
                << token::TAB << rmsOPD
                << token::TAB << rmsOPDtilt
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //
