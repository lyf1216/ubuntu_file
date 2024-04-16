/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    setWaveField

Description
    Loop over every cell in the computational domain and set VOF-ratio and
    velocity field accordingly to specified wave theory.

Author
    Niels Gjoel Jacobsen, Technical University of Denmark.  All rights reserved.

Additional information
    Implementation published and validated in the following journal article:

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "fvMesh.H"
#include "volFields.H"
#include "setWaveField.H"

#include "uniformDimensionedFields.H"

// #include "crossVersionCompatibility.H"
// #include "externalWaveForcing.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

#   include "readGravitationalAcceleration.H"

    Info << "\nReading waveProperties" << endl;

    IOdictionary waveProperties
    (
        IOobject
        (
            "waveProperties_OF+",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    /* Compute a reference point that is placed on the free surface */
    // A zero offset between original and new reference level. Default offset
    dimensionedVector referencePoint("NULL", dimLength, vector::zero);

	// Get the water level
    scalar waterLevel = 0;
    if (waveProperties.found("INLET"))
    {
	    waterLevel = waveProperties.subDict("INLET").get<scalar>("waterDepthRef");
    }
    else if (waveProperties.found("inlet"))
    {
	    waterLevel = waveProperties.subDict("inlet").get<scalar>("waterDepthRef");
    }

    referencePoint.value() = g.value()/Foam::mag(g.value());
    referencePoint.value() = Foam::cmptMag(referencePoint.value());
    referencePoint.value() *= waterLevel;


    Info<< "\nReading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            Foam::waves2Foam::aName(),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            Foam::waves2Foam::pName(),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Setting the wave field ...\n" << endl;
    IOdictionary transProp
    (
        IOobject
        (
            "transportProperties_dim",
            "constant",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& sDA(transProp.subDict(Foam::waves2Foam::airPhase()));
    const dimensionedScalar rho2("rho_air", dimDensity, sDA.get<scalar>("rho"));
    
    const dictionary& sDW(transProp.subDict(Foam::waves2Foam::waterPhase()));
    const dimensionedScalar rho1("rho_water", dimDensity, sDW.get<scalar>("rho"));

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha*rho1 + (1-alpha)*rho2
    );

    //----------------结果输出
 
    volScalarField gh("gh", g & (mesh.C() - referencePoint));

    Info<< "Reading field p\n" << endl;

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd + rho*gh
    );

    // rho.write();
    p.write();

    Info << nl << "End" << nl << endl;
    return 0;
}