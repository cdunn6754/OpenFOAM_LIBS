/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::functionObjects::customCloudInfo

Description
    This function object outputs Lagrangian cloud information to a file. It is best
    with uncoupledCoalChemistryFoam or similar (e.g. uncoupledKinematicParcelFoam). It 
    tracks all of the parcels of a particular cloud (hard-coded coalCloud1). Then the 
    parcel variables of interest (e.g. Temp, Ygas) are written to a file for every
    timestep. The current varibles available are T, Ygas, Ysolid, mass0. When writing 
    the names of the variables in the dictionary write them as they appear in the 
    normal solver output files.

    Example of function object specification, put this in the controlDict:
    \verbatim
    functions
    {
    cloudInfo1
        {
          type customCloudInfo;
          libs ("libcustomLagrangianFunctionObjects.so");
          parcelProperties
	    (
	     T
	     Ysolid
	     Ygas
	     );
	 }
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::writeFiles

SourceFiles
    customCloudInfo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_customCloudInfo_H
#define functionObjects_customCloudInfo_H

#include "writeFiles.H"

// Clint 06-14-17
#include "kinematicCloud.H"
#include "basicReactingMultiphaseCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class customCloudInfo Declaration
\*---------------------------------------------------------------------------*/

class customCloudInfo
:
    public writeFiles
{
protected:

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        customCloudInfo(const customCloudInfo&);

        //- Disallow default bitwise assignment
        void operator=(const customCloudInfo&);

  // Clint 06-14-17
  // Adding pointer to the cloud as a member so we can access info in all functions
  basicReactingMultiphaseCloud* cloud;

public:

    //- Runtime type information
    TypeName("customCloudInfo");


    // Constructors

        //- Construct from Time and dictionary
        customCloudInfo
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~customCloudInfo();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
