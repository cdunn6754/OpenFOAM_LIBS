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

\*---------------------------------------------------------------------------*/

#include "parcelPropertyInfo.H"
#include "addToRunTimeSelectionTable.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(parcelPropertyInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        parcelPropertyInfo,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::parcelPropertyInfo::writeFileHeader(const label i)
{
  // need to know how many parcels so we can write the headers for them
  label nParcels = cloud->nParcels();

  writeHeader(files()[i], "Parcel " + names()[i]);
  files()[i] <<  "Time" << token::TAB << token::TAB ;
      // loop through all parcels and make their column headers
      for (int j=0; j < nParcels; j++)
	{
	  int parcelNumber = j;
	  files()[i]
	    << token::TAB << "p" << parcelNumber << token::TAB << token::TAB;
	}
      files()[i] << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::parcelPropertyInfo::parcelPropertyInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
  writeFiles(name, runTime, dict, name)
{
    read(dict);
    //const word& cloudName = names()[i];
    const word& cloudName = "coalCloud1";

    const kinematicCloud& cloud1 =
      obr_.lookupObject<kinematicCloud>(cloudName);

    // Type-cast to something we can work with
    cloud = (basicReactingMultiphaseCloud*) &cloud1;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::parcelPropertyInfo::~parcelPropertyInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::parcelPropertyInfo::read(const dictionary& dict)
{
    writeFiles::resetNames(dict.lookup("parcelProperties"));

    Info<< type() << " " << name() << ": ";
    if (names().size())
    {
        Info<< "applying to lagrangian fields:" << nl;
        forAll(names(), i)
        {
            Info<< "    " << names()[i] << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no lagrangian fields to be processed" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::parcelPropertyInfo::execute()
{
    return true;
}


bool Foam::functionObjects::parcelPropertyInfo::write()
{
  writeFiles::write();
      // total number of parcels
      label parcelCount = cloud->first()->particleCount_;

      // iterate through requested parcel properties (e.g. T, mass ...)
      forAll(names(), i)
	{
	  // scalarField to hold the parcel property
	  // one entry per parcel (even if they are deleted from cloud)
	  scalarField parcelFieldValues(parcelCount,0.);

	  // paired to parcelFieldValues, carries information about 
	  // the parcels existence in the cloud
	  scalarField parcelExistence(parcelCount,0.);

	  // Stuff used for pulling composition info from parcels
	  label idGas = cloud->composition().idGas();
	  label idSolid = cloud->composition().idSolid();

	  // Iterate through all parcels still in the cloud
	  // any parcels already deleted from the cloud
	  // will have their values remain at -1.0
	  forAllIter(basicReactingMultiphaseCloud, *cloud,  iter)
	    { 
	      // original parcel ID, this will serve as an iterator 
	      // through parcelFieldValues, will be unique as long
	      // as it's run in parallel (or at least the parcels
	      // are all  injected into the same processor)
	      label parcelOrigId = iter().origId();


	      // Could not figure out how to reflect name string to function name
	      // Now we use the dumb way to get from string "fieldName" i.e.
	      // "T" to the access function iter().T()
	      if (names()[i] == "T")
	      	{
	      	  parcelFieldValues[parcelOrigId] = 
	      	    returnReduce(iter().T(), sumOp<scalar>());
	      	}

	      else if (names()[i] == "Ygas")
	      	{
	      	  parcelFieldValues[parcelOrigId] = 
	      	    returnReduce(iter().Y()[idGas], sumOp<scalar>());
	      	}

	      else if (names()[i] == "Ysolid")
	      	{
	      	  parcelFieldValues[parcelOrigId] = 
	      	    returnReduce(iter().Y()[idSolid], sumOp<scalar>());
	      	}

	      else if (names()[i] == "mass0")
	      	{
	      	  parcelFieldValues[parcelOrigId] = 
	      	    returnReduce(iter().mass0(), sumOp<scalar>());
	      	}

	      else if (names()[i] == "mass")
	      	{
	      	  parcelFieldValues[parcelOrigId] = 
	      	    returnReduce(iter().mass(), sumOp<scalar>());
	      	}

	      else if (names()[i] == "Yash")
	      	{
	      	  parcelFieldValues[parcelOrigId] = 
	      	    returnReduce(iter().YSolid()[1]*iter().Y()[idSolid], 
				 sumOp<scalar>());
	      	}

	      else if (names()[i] == "YC")
	      	{
	      	  parcelFieldValues[parcelOrigId] = 
	      	    returnReduce(iter().YSolid()[0]*iter().Y()[idSolid],
				 sumOp<scalar>());
	      	}

	      Info << "Solid fraction: " << iter().Y()[idSolid] << endl; 
	      
	    }


	  if (Pstream::master())
	    {
	      writeTime(files()[i]);
	      files()[i]
		<< token::TAB;
	      // write out the fields for every parcel
	      for (label k = 0; k < parcelCount; k++)
		{
		  if (parcelExistence[k] == 1)
		    {
		      files()[i]
			<< 
			parcelFieldValues[k]
			<< token::TAB;
		    }
		  else // if the parcel is deleted, write NaN
		    {
		      files()[i]
			<< 
			"NaN"
			<< token::TAB;
		    }
		  
		}
	      files()[i]
		<< endl;
	    }
	}

  return true;
}


// ************************************************************************* //
