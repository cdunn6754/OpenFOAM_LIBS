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
	  int parcelNumber = j + 1;
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

  forAll(names(), i)
    {
      // need to know how many parcels
      label nParcels = cloud->nParcels();
      
      scalarField parcelFieldValues(nParcels,0.);

      label idGas = cloud->composition().idGas();
      //label idLiquid = cloud->composition().idLiquid();
      label idSolid = cloud->composition().idSolid();



      label j = 0;
      forAllIter(basicReactingMultiphaseCloud, *cloud,  iter)
	{

	  if (iter().active())
	    {
	      // Could not figure out how to reflect name string to function name
	      // Now we use the dumb way to get from string "fieldName" i.e.
	      // "T" to the access function iter().T()
	      if (names()[i] == "T")
		{
		  parcelFieldValues[j] = 
		    returnReduce(iter().T(), sumOp<scalar>());
		}

	      else if (names()[i] == "Ygas")
		{
		  parcelFieldValues[j] = 
		    returnReduce(iter().Y()[idGas], sumOp<scalar>());
		}

	      else if (names()[i] == "Ysolid")
		{
		  parcelFieldValues[j] = 
		    returnReduce(iter().Y()[idSolid], sumOp<scalar>());
		}

	      else if (names()[i] == "mass0")
		{
		  parcelFieldValues[j] = 
		    returnReduce(iter().mass0(), sumOp<scalar>());
		}

	      else if (names()[i] == "mass")
		{
		  parcelFieldValues[j] = 
		    returnReduce(iter().mass(), sumOp<scalar>());
		}

	      j++;
	    }

	  else // if the parcel is no longer active
	    {
	      parcelFieldValues[j] = 0.0;
	    }
	  // double check that we don't overrun the cloud
	  if (j>= nParcels)
	    {
	      break;
	    }
	      
	}


      if (Pstream::master())
        {
	  writeTime(files()[i]);
	  files()[i]
	    << token::TAB;
	  // write out the fields for every parcel
	  for (label k = 0; k < nParcels; k++)
	    {
	      files()[i]
		<< parcelFieldValues[k] << token::TAB;
	    }
	  files()[i]
	    << endl;
	}
    }

  return true;
}


// ************************************************************************* //
