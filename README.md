# laserOPD

Particle-based library for calculating the optical path difference (OPD) of a laser beam through a non-uniform density field via OpenFOAM.

Detailed information about the theoretical background and the implementation can 
be found in:

 [Manuscript for softwareX (laserOPD: particle-based library for calculating the optical path difference of a laser beam via OpenFOAM)](documentation\main.pdf) 

## Supported OpenFOAM Versions:

 * OpenFOAM (ORG) v7 - v8
 * OpenFOAM (ESI) v2006 - v2106

## Authors

 * Xiang Ren <renxiangbuaa@gmail.com>

## Installation

1. Clone branch corresponding to the version of OpenFOAM used

   ```shell
   $ git clone https://github.com/renxiangbuaa/laserOPD.git
   ```

2. Build the `laserOPD` and `laserOPDInfoFunctionObjects` library

   ```shell
   $ ./Allwmake
   ```



## Usage

To use the `laserOPD` library you have to add the library to your controlDict by editing `system/controlDict`

```c++
libs("liblaserOPD.so");
```

Within adding `constant/radiationProperties` file

```c++
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7  
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      radiationProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

radiation       on;

// Selection of laserOPD models
radiationModel  laserOPD;

// Number of flow iterations per radiation iteration
solverFreq 1;

// Model for laser Beam
laserModel             circularLaser;//squareLaser;//circularInSqureLaser
  circularLaserCoeffs
  {
    nTheta	              100;
    nr 	                  100;
    radius                0.03;
  }

  squareLaserCoeffs
  {
    nLength               100;
    length                0.08;
    directionX            (1 0 0);
  }

  circularInSquareLaserCoeffs
  {
    nLength               100;
    length                0.08;
    directionX            (1 0 0);
    radius                0.03;
  }

  // Maximum tracking length of the laser beam
  maxTrackLength        0.5;
  
  // Centre of the laser beam
  focalLaserPosition    (0 0.5 0);

  // Direction of the laser beam
  laserDirection        (0 -1 0);

  // Gladstone-Dale constant
  // which is based on gas composition and laser wavelength
  GladstoneDale         0.00022892515122;

  // Reference density
  densityRef            0.9575
```

Further, add the following to `system/controlDict` to count the average optical path length (OPL), the tilt factor and the root-mean-square of the optical path difference (OPD)

```c++
functions
{
    OPDInfo
    {
        type            laserOPDInfo;
        libs            (laserOPDInfoFunctionObjects); // for of-ESI v2006-v2106
        // libs            ("liblaserOPDInfoFunctionObjects.so"); // for of-ORG v7-v8
        writeControl    timeStep;
        writeInterval   1;
        clouds
        (
            laserCloud
        );
    }
}
```



## Tutorial

The `laserOPD` library is inherited from `radiationModels` in OpenFOAM, so the radiation model needs to be added to the corresponding solver first. 

Here we use the compressible solver `rhoPimpleFOAM` as an example to explain how to add `laserOPD` . 

1. Add three codes to `rhoPimpleFOAM.C`,

```c++
...
// addForlaserOPD: radiation model header file
#include "radiationModel.H"
...
  
int main(int argc, char *argv[])
{
    ...
    #include "createFields.H"

    // addForlaserOPD: Creating radiation model
    // for of-v2006 and of-v2106
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(rho)
    );
    // // for of-v7 and of-v8
    // autoPtr<radiationModel> radiation
    // (
    //     radiationModel::New(rho)
    // );
    ...

    while (runTime.run())
    {
        ...
          
        rho = thermo.rho();
        
        // addForlaserOPD: loop
        radiation->correct();
        ...
    }
  ...
}
```

2. add the `radiationModels` library to `Make/options`,

The modified solver code can be found in the folder `myRhoPimpleFoam`, which can be compiled to produce a new solver by running `wmake`.



## Tests

A flow over cylinder case with Re=3900 and Ma=0.4 can be used to test the new solver `myRhoPimpleFoam` with `laserOPD` library.




## License 

This OpenFOAM library is under the GNU General Public License. 
