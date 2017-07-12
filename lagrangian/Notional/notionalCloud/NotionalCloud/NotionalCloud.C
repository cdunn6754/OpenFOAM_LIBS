/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "NotionalCloud.H"

#include "DevolatilisationModel.H"
#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NotionalCloud<CloudType>::NotionalCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    notionalCloud(),
    cloudCopyPtr_(NULL)
{
    if (this->solution().active())
    {
        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
      //resetSourceTerms();
    }
}


template<class CloudType>
Foam::NotionalCloud<CloudType>::NotionalCloud
(
    NotionalCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    notionalCloud(),
    cloudCopyPtr_(NULL)
{}


template<class CloudType>
Foam::NotionalCloud<CloudType>::NotionalCloud
(
    const fvMesh& mesh,
    const word& name,
    const NotionalCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    notionalCloud(),
    cloudCopyPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NotionalCloud<CloudType>::~NotionalCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




template<class CloudType>
void Foam::NotionalCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<NotionalCloud<CloudType>> td(*this);

        this->solve(td);
    }
}


template<class CloudType>
void Foam::NotionalCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    typedef typename particle::TrackingData<NotionalCloud<CloudType>>
        tdType;

    tdType td(*this);

    Cloud<parcelType>::template autoMap<tdType>(td, mapper);

    this->updateMesh();
}



template<class CloudType>
void Foam::NotionalCloud<CloudType>::writeFields() const
{
    if (this->size())
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


// ************************************************************************* //
