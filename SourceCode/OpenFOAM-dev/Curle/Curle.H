/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::Curle

Description
    Implementation of the Curle analogy for incompressible flow as
    described in:
    \verbatim
        Johan Larsson, Lars Davidson, Magnus Olsson, Lars-Erik Eriksson (2004).
        Aeroacoustic Investigation of an Open Cavity at Low Mach Number.
        AIAA Journal, Vol. 42, No. 12.
    \endverbatim

    Example of function object specification:
    \verbatim
    Curle1
    {
        type                    Curle;
        functionObjectLibs      ("libAcousticAnalogy.so");
        log                     true;

        patches                 ( cylinder1 cylinder2 );
        cellZone                cellZoneName;
        rhoRef                  1.225;
        cRef                    343;

        observers
        {
            microphone-A
            {
                position    	(-0.4760595 1.58962725 0);
            }
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description                                       | Required    | Default value
        type         | type name: Curle                                  | yes         |
        log          | write acoustic data to standard output            | no          | no
        patches      | patches included in the acoustic calculation      | yes         |
        cellZone     | name of cell zone for acoustic volume integrals   | no          |
        rhoRef       | reference density                                 | yes         |
        cRef         | reference speed of sound                          | yes         |
        pName        | pressure field name                               | no          | p
        UName        | velocity field name                               | no          | U
    \endtable


\*---------------------------------------------------------------------------*/

#ifndef Curle_H
#define Curle_H

#include "regionFunctionObject.H"
#include "logFiles.H"
#include "volFields.H"
#include "fvcDdt.H"
#include "soundObserver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class Curle Declaration
\*---------------------------------------------------------------------------*/

class Curle
:
    public regionFunctionObject,
    public logFiles
{
protected:

    // Protected data

        //- Name of this set of Curle,
        word name_;

        //- Initialised flag
        bool initialised_;
        
        //- Switch to send output to Info as well as to file
        Switch log_;
        
        //- Patches to sample
        labelHashSet patches_;
        
        //- Name of cellZone
        word cellZoneName_;
        
        //- ID of cellZone
        label cellZoneID_;
        
        //- Name of pressure field (optional)
        word pName_;
        
        //- Name of velocity field (optional)
        word UName_;

        //- Reference density needed for incompressible cases
        scalar rhoRef_;
        
        //- Reference speed of sound
        scalar cRef_;

        //- List of observers
        List<SoundObserver> observers_;


    // Protected Member Functions

        //- Create file names for Curle
        wordList createFileNames(const dictionary& dict) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Helper function to write acoustic data
        void writeCurle();

        //- Initialise the fields
        void initialise();

        //- Return pressure
        tmp<volScalarField> p() const;

        //- Return pressure time derivative
        tmp<volScalarField> dpdt() const;

        //- Return Lighthill tensor
        tmp<volTensorField> Tij() const;

        //- Return first time derivative of Lighthill tensor 
        tmp<volTensorField> dTijdt() const;

        //- Return second time derivative of Lighthill tensor 
        tmp<volTensorField> d2Tijdt2() const;

        //- Disallow default bitwise copy construct
        Curle(const Curle&);

        //- Disallow default bitwise assignment
        void operator=(const Curle&);


public:

    //- Runtime type information
    TypeName("Curle");


    // Constructors

        //- Construct from Time and dictionary
        Curle
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
        
        //- Construct from objectRegistry and dictionary
        Curle
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


    //- Destructor
    virtual ~Curle();


    // Member Functions

        //- Read the Curle data
        virtual bool read(const dictionary&);

        //- Calculate the acoustic data
        virtual void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the acoustic data
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
