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

#include "Curle.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Curle, 0);
    addToRunTimeSelectionTable(functionObject, Curle, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::Curle::createFileNames(const dictionary& dict) const
{
    DynamicList<word> names(1);

    const word CurleType(dict.lookup("type"));

    names.append(CurleType);

    return names;
}


void Foam::functionObjects::Curle::writeFileHeader(const label i)
{
    if (i == 0)
    {
        writeHeader(file(i), "Curle Acoustic Analogy");
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


void Foam::functionObjects::Curle::writeCurle()
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


void Foam::functionObjects::Curle::initialise()
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
    
    if
    (
        cellZoneID_ == -1
     && cellZoneName_ != word::null    
    )
    {
        FatalErrorInFunction
            << "Could not find cellZone '" << cellZoneName_ << "' in database." << nl
            << exit(FatalError);
    }

    initialised_ = true;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::Curle::p() const
{
    return
    (
        rhoRef_ * obr_.lookupObject<volScalarField>(pName_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::Curle::dpdt() const
{
    return
    (
        rhoRef_ * Foam::fvc::ddt(obr_.lookupObject<volScalarField>(pName_))
    );
}


Foam::tmp<Foam::volTensorField> Foam::functionObjects::Curle::Tij() const
{
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

    return
    (
        rhoRef_*(U*U)
    );
}


Foam::tmp<Foam::volTensorField> Foam::functionObjects::Curle::dTijdt() const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Get velocity field of current and last two time steps
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volVectorField& U0 = obr_.lookupObject<volVectorField>(UName_).oldTime();
    const volVectorField& U00 = obr_.lookupObject<volVectorField>(UName_).oldTime().oldTime();

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
        rhoRef_*
        (
            (
                coefft*U*U
              - coefft0*U0*U0
              + coefft00*U00*U00
            )
          / mesh.time().deltaT()
        )
    );
}


Foam::tmp<Foam::volTensorField> Foam::functionObjects::Curle::d2Tijdt2() const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Get velocity field of current and last two time steps
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volVectorField& U0 = obr_.lookupObject<volVectorField>(UName_).oldTime();
    const volVectorField& U00 = obr_.lookupObject<volVectorField>(UName_).oldTime().oldTime();

    // Get time stepping size
    scalar deltaT = mesh.time().deltaTValue();
    scalar deltaT0 = mesh.time().deltaT0Value();
    scalar dt = 4.0 / sqr(deltaT + deltaT0);

    // Calculate coefficients
    scalar coefft   = (deltaT + deltaT0)/(2.*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2.*deltaT0);
    scalar coefft0  = coefft + coefft00;

    // First-order Euler second time derivative
    // see: EulerD2dt2Scheme.H
    return
    (
        rhoRef_*dt*
        (
            (
                coefft*U*U
              - coefft0*U0*U0
              + coefft00*U00*U00
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::Curle::Curle
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    name_(name),
    initialised_(false),
    log_(false),
    patches_(0),
    cellZoneName_(word::null),
    cellZoneID_(-1),
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


Foam::functionObjects::Curle::Curle
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    regionFunctionObject(name, obr, dict),
    logFiles(obr_, name),
    name_(name),
    initialised_(false),
    log_(false),
    patches_(0),
    cellZoneName_(word::null),
    cellZoneID_(-1),
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


Foam::functionObjects::Curle::~Curle()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::Curle::read(const dictionary& dict)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    log_ = dict.lookupOrDefault<Switch>("log", false);

    patches_ = mesh.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));
    
    cellZoneName_ = dict.lookupOrDefault<word>("cellZone", word::null);
    
    cellZoneID_ = mesh.cellZones().findZoneID(cellZoneName_);
    
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


bool Foam::functionObjects::Curle::execute()
{
    return true;
}


bool Foam::functionObjects::Curle::write()
{
    calculate();

    if (Pstream::master())
    {
        logFiles::write();

        writeCurle();

        if (log_) Info << endl;
    }
    
    return true;
}


void Foam::functionObjects::Curle::calculate()
{
    initialise();

    // Get access to mesh
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    
    // Pressure field and its time derivative
    volScalarField p = this->p();
    p.correctBoundaryConditions();
    volScalarField dpdt = this->dpdt();
    dpdt.correctBoundaryConditions();

    // Lighthill tensor and its time derivatives
    volTensorField Tij = this->Tij();
    volTensorField dTijdt = this->dTijdt();
    volTensorField d2Tijdt2 = this->d2Tijdt2();

    // Calculate constant
    scalar coeff = 1.0 / ( 4.0 * Foam::constant::mathematical::pi );

    // Identity tensor
    tensor I(1,0,0,0,1,0,0,0,1);

    // Loop over all observer
    forAll(observers_, obsI)
    {
        SoundObserver& obs = observers_[obsI];
        scalar pPrime = 0.0;

        // Volume integral
        if (cellZoneID_ != -1)
        {
            // List of cells in cellZoneID
            const labelList& cells = mesh.cellZones()[cellZoneID_];

            // Cell volume and cell center
            const scalarField& V = mesh.V();
            const vectorField& C = mesh.C();

            // Loop over all cells
            forAll(cells, i)
            {
                label cellI = cells[i];

                // Distance to observer
                scalar r = mag(obs.position() - C[cellI]);
                vector l = (obs.position() - C[cellI]) / r;

                // Calculate pressure fluctuation
                pPrime += coeff * V[cellI] * 
                        (
                            ((l*l) && d2Tijdt2[cellI]) / (cRef_ * cRef_ * r)
                          + ((3.0 * l*l - I) && dTijdt[cellI]) / (cRef_ * r * r)
                          + ((3.0 * l*l - I) && Tij[cellI]) / (r * r * r)
                        );
            }
            reduce(pPrime, sumOp<scalar>());
        }



        // Surface integral - loop over all patches
        forAllConstIter(labelHashSet, patches_, iter)
        {
        	// Get patch ID
            label patchI = iter.key();

            // Surface area vector and face center at patch
            vectorField Sf = mesh.Sf().boundaryField()[patchI];
            vectorField Cf = mesh.Cf().boundaryField()[patchI];
            
            // Normal vector pointing towards fluid
            vectorField n = -Sf/mag(Sf);

            // Pressure field and time derivative at patch
            scalarField pp = p.boundaryField()[patchI];
            scalarField dpdtp = dpdt.boundaryField()[patchI];
            
            // Lighthill tensor on patch
            tensorField Tijp = Tij.boundaryField()[patchI];
            tensorField dTijdtp = dTijdt.boundaryField()[patchI];

            // Distance surface-observer
            scalarField r = mag(obs.position() - Cf);
            vectorField l = (obs.position() - Cf) / r;

            // Calculate pressure fluctuations
            pPrime += coeff * gSum
            (
                (
                    (l*n)
                  && 
                    (
                        (dpdtp*I - dTijdtp) / (cRef_*r)
                      + (pp*I - Tijp) / sqr(r)
                    )
                )
              * mag(Sf)
            );
        }
        obs.pPrime(pPrime);
    }
}


// ************************************************************************* //
