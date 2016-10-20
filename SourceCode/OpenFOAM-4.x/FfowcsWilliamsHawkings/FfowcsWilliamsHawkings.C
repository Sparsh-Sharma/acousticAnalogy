/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "FfowcsWilliamsHawkings.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(FfowcsWilliamsHawkings, 0);
    addToRunTimeSelectionTable(functionObject, FfowcsWilliamsHawkings, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::FfowcsWilliamsHawkings::createFileNames(const dictionary& dict) const
{
    DynamicList<word> names(1);

    const word FfowcsWilliamsHawkingsType(dict.lookup("type"));

    names.append(FfowcsWilliamsHawkingsType);

    return names;
}


void Foam::functionObjects::FfowcsWilliamsHawkings::writeFileHeader(const label i)
{
    if (i == 0)
    {
        writeHeader(file(i), "FfowcsWilliamsHawkings Acoustic Analogy");
        writeCommented(file(i), "Time");

        forAll(observers_, obsI)
        {
            file(i)
	        << "pPrime at " << observers_[obsI].name() << "   ";
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unhandled file index: " << i
            << abort(FatalError);
    }

    file(i)<< endl;
}


void Foam::functionObjects::FfowcsWilliamsHawkings::writeFfowcsWilliamsHawkings()
{
    // Log output
    if (log_) Info
        << type() << " output:" << endl;
        
    forAll(observers_, obsI)
    {
        if (log_) Info
            << observers_[obsI].name() << "   " << observers_[obsI].pPrime() << endl;
    }

    // File output
    file(0) << obr_.time().value() << tab << setw(1) << "   ";
    forAll(observers_, obsI)
    {
        file(0)
            << observers_[obsI].pPrime() << "   ";
    }
    file(0) << endl;
}


void Foam::functionObjects::FfowcsWilliamsHawkings::initialise()
{
    if (initialised_)
    {
        return;
    }

    if
    (
        !obr_.foundObject<volVectorField>(UName_)
     || !obr_.foundObject<volScalarField>(pName_)
    )
    {
        FatalErrorInFunction
            << "Could not find " << UName_ << ", " << pName_
            << exit(FatalError);
    }

    initialised_ = true;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::FfowcsWilliamsHawkings::p() const
{
    return
    (
        rhoRef_ * obr_.lookupObject<volScalarField>(pName_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::FfowcsWilliamsHawkings::dpdt() const
{
    return
    (
        rhoRef_ * Foam::fvc::ddt(obr_.lookupObject<volScalarField>(pName_))
    );
}


Foam::tmp<Foam::volVectorField> Foam::functionObjects::FfowcsWilliamsHawkings::U() const
{
    return
    (
        obr_.lookupObject<volVectorField>(UName_)
    );
}


Foam::tmp<Foam::volVectorField> Foam::functionObjects::FfowcsWilliamsHawkings::dUdt() const
{
    return
    (
        Foam::fvc::ddt(obr_.lookupObject<volVectorField>(UName_))
    );
}


Foam::tmp<Foam::surfaceVectorField> Foam::functionObjects::FfowcsWilliamsHawkings::n() const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    return
    (
        -mesh.Sf() / mesh.magSf()
    );
}


Foam::tmp<Foam::surfaceVectorField> Foam::functionObjects::FfowcsWilliamsHawkings::dndt() const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Get surface normal vector of current and last two time steps
    surfaceVectorField n = -mesh.Sf() / mesh.magSf();
    surfaceVectorField n0 = -mesh.Sf().oldTime() / mesh.magSf().oldTime();
    surfaceVectorField n00 = -mesh.Sf().oldTime().oldTime() / mesh.magSf().oldTime().oldTime();

    // Get time stepping size
    scalar deltaT = mesh.time().deltaTValue();
    scalar deltaT0 = mesh.time().deltaT0Value();

    // Calculate coefficients
    scalar coefft   = 1.0 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    // Second-order backward-differencing time derivative
    // assuming constant time stepping
    // see: backwardDdtScheme.C
    return
    (
        (
            (
                coefft*n
              - coefft0*n0
              + coefft00*n00
            )
          / mesh.time().deltaT()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::FfowcsWilliamsHawkings::FfowcsWilliamsHawkings
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    writeFiles(name, runTime, dict, name),
    name_(name),
    initialised_(false),
    log_(false),
    patches_(0),
    pName_(word::null),
    UName_(word::null),
    rhoRef_(-1),
    cRef_(-1),
    observers_(0)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
    resetNames(createFileNames(dict));
}


Foam::functionObjects::FfowcsWilliamsHawkings::FfowcsWilliamsHawkings
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    writeFiles(name, obr, dict, name),
    name_(name),
    initialised_(false),
    log_(false),
    patches_(0),
    pName_(word::null),
    UName_(word::null),
    rhoRef_(-1),
    cRef_(-1),
    observers_(0)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
    resetNames(createFileNames(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::functionObjects::FfowcsWilliamsHawkings::~FfowcsWilliamsHawkings()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::FfowcsWilliamsHawkings::read(const dictionary& dict)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    log_ = dict.lookupOrDefault<Switch>("log", false);

    patches_ = mesh.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));
    
    pName_ = dict.lookupOrDefault<word>("pName", "p");
    
    UName_ = dict.lookupOrDefault<word>("UName", "U");
    
    rhoRef_ = readScalar(dict.lookup("rhoRef"));
    
    cRef_ = readScalar(dict.lookup("cRef"));

    if (log_) Info << "Acoustic analogy: " << type() << nl << nl;
        
    // Read observers
    {
        const dictionary& obsDict = dict.subDict("observers");
        wordList obsNames = obsDict.toc();

        forAll (obsNames, obsI)
        {
            word oName = obsNames[obsI];
            vector oPos (vector::zero);
            obsDict.subDict(oName).lookup("position") >> oPos;

            observers_.append
            (
                SoundObserver
                (
                    oName,
                    oPos
                )
            );
        }
    }
    
    return true;
}


bool Foam::functionObjects::FfowcsWilliamsHawkings::execute()
{
    return true;
}


bool Foam::functionObjects::FfowcsWilliamsHawkings::write()
{
    calculate();

    if (Pstream::master())
    {
        writeFiles::write();

        writeFfowcsWilliamsHawkings();

        if (log_) Info << endl;
    }
    
    return true;
}


void Foam::functionObjects::FfowcsWilliamsHawkings::calculate()
{
    initialise();

    // Get access to mesh
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    
    // Pressure field and its time derivative
    volScalarField p = this->p();
    volScalarField dpdt = this->dpdt();
    
    // Velocity field and its time derivative
    volVectorField U = this->U();
    volVectorField dUdt = this->dUdt();

    // Surface normal vector and its time derivative
    surfaceVectorField n = this->n();
    surfaceVectorField dndt = this->dndt();

    // Calculate constant
    scalar coeff = 1.0 / ( 4.0 * Foam::constant::mathematical::pi );

    // Identity tensor
    tensor I(1,0,0,0,1,0,0,0,1);

    // Loop over all observer
    forAll(observers_, obsI)
    {
        SoundObserver& obs = observers_[obsI];
        scalar pPrime = 0.0;

        // Surface integral - loop over all patches
        forAllConstIter(labelHashSet, patches_, iter)
        {
        	// Get patch ID
            label patchI = iter.key();

            // Surface area vector and face center at patch
            scalarField magSf = mesh.magSf().boundaryField()[patchI];
            vectorField Cf = mesh.Cf().boundaryField()[patchI];
            
            // Normal vector pointing towards fluid
            vectorField np = n.boundaryField()[patchI];
            vectorField dndtp = dndt.boundaryField()[patchI];

            // Pressure field and time derivative at patch
            scalarField pp = p.boundaryField()[patchI];
            scalarField dpdtp = dpdt.boundaryField()[patchI];
            
            // Velocity field ant time derivative at patch
            vectorField Up = U.boundaryField()[patchI];
            vectorField dUdtp = dUdt.boundaryField()[patchI];

            // Distance surface-observer
            scalarField r = mag(obs.position() - Cf);
            vectorField l = (obs.position() - Cf) / r;
            vectorField rHat = l / r;
            
            // Surface normal velocity
            scalarField Un = Up & np;
            scalarField UnDot = (dUdtp & np) + (Up & dndtp);
            
            // Mach number
            vectorField M = Up / cRef_;
            scalarField Mr = M & rHat;
            scalarField Mn = M & np;
            scalarField MrDot = rHat & (dUdtp / cRef_);

            // Intermediate variables
            vectorField Li = (p.boundaryField()[patchI]*I) & np;
            vectorField dLidt = (dpdt.boundaryField()[patchI]*I) & np;

            // Thickness noise
            pPrime += gSum( UnDot * magSf / (r * sqr(1.0 - Mr)) ) * rhoRef_ * coeff;
            pPrime += gSum( Un * (r*MrDot + cRef_*(Mr - magSqr(M))) / ( sqr(r) * sqr(1.0 - Mr) * (1.0 - Mr) ) ) * rhoRef_ * coeff;
            
            // Loading noise
            pPrime += gSum( (dLidt & rHat) * magSf / (r * sqr(1.0 - Mr)) ) * coeff / cRef_;
            pPrime += gSum( ( (Li & rHat) - (Li & M) ) * magSf / (sqr(r) * sqr(1.0 - Mr)) ) * coeff;
            pPrime += gSum( ( (Li & rHat) * (r * MrDot + cRef_*(Mr - magSqr(M))) * magSf) / (sqr(r) * sqr(1.0 - Mr) * (1.0 - Mr)) ) * coeff;
        }
        obs.pPrime(pPrime);
    }
}


// ************************************************************************* //
