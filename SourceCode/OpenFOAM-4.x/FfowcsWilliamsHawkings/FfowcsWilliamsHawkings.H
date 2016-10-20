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
    Foam::FfowcsWilliamsHawkings

Description
    Implementation of the Ffowcs-Williams-Hawkings analogy for incompressible flow as
    described in:
    \verbatim
        F. Farassat (2007).
        Derivation of Formulations of 1 and 1A of Farassat.
        NASA Technical Memorandum 214853.
    \endverbatim

    Example of function object specification:
    \verbatim
    FfowcsWilliamsHawkings1
    {
        type                    FfowcsWilliamsHawkings;
        functionObjectLibs      ("libAcousticAnalogy.so");
        log                     true;

        patches                 ( cylinder1 cylinder2 );
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
        type         | type name: FfowcsWilliamsHawkings                                    | yes         |
        log          | write acoustic data to standard output            | no          | no
        patches      | patches included in the acoustic calculation      | yes         |
        rhoRef       | reference density                                 | yes         |
        cRef         | reference speed of sound                          | yes         |
        pName        | pressure field name                               | no          | p
        UName        | velocity field name                               | no          | U
    \endtable


\*---------------------------------------------------------------------------*/

#ifndef FfowcsWilliamsHawkings_H
#define FfowcsWilliamsHawkings_H

#include "writeFiles.H"
#include "volFields.H"
#include "fvcDdt.H"
#include "soundObserver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class FfowcsWilliamsHawkings Declaration
\*---------------------------------------------------------------------------*/

class FfowcsWilliamsHawkings
:
    public writeFiles
{
protected:

    // Protected data

        //- Name of this set of FfowcsWilliamsHawkings,
        //  Also used as the name of the probes directory.
        word name_;

        //- Initialised flag
        bool initialised_;
        
        //- Switch to send output to Info as well as to file
        Switch log_;
        
        //- Patches to sample
        labelHashSet patches_;
        
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

        //- Create file names for FfowcsWilliamsHawkings
        wordList createFileNames(const dictionary& dict) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Helper function to write acoustic data
        void writeFfowcsWilliamsHawkings();

        //- Initialise the fields
        void initialise();

        //- Return pressure field
        tmp<volScalarField> p() const;

        //- Return pressure field time derivative
        tmp<volScalarField> dpdt() const;

        //- Return velocity field
        tmp<volVectorField> U() const;

        //- Return velocity field time derivative
        tmp<volVectorField> dUdt() const;

        //- Return surface normal vector
        tmp<surfaceVectorField> n() const;

        //- Return surface normal vector time derivative
        tmp<surfaceVectorField> dndt() const;

        //- Disallow default bitwise copy construct
        FfowcsWilliamsHawkings(const FfowcsWilliamsHawkings&);

        //- Disallow default bitwise assignment
        void operator=(const FfowcsWilliamsHawkings&);


public:

    //- Runtime type information
    TypeName("FfowcsWilliamsHawkings");


    // Constructors

        //- Construct from Time and dictionary
        FfowcsWilliamsHawkings
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
        
        //- Construct from objectRegistry and dictionary
        FfowcsWilliamsHawkings
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


    //- Destructor
    virtual ~FfowcsWilliamsHawkings();


    // Member Functions

        //- Read the FfowcsWilliamsHawkings data
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
