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

#include "DaemReactingMultiphaseCloud.H"

#include "DaemDevolatilisationModel.H"
#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::setModels()
{
    devolatilisationModel_.reset
    (
        DaemDevolatilisationModel<DaemReactingMultiphaseCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    surfaceReactionModel_.reset
    (
        SurfaceReactionModel<DaemReactingMultiphaseCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::cloudReset
(
    DaemReactingMultiphaseCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    devolatilisationModel_.reset(c.devolatilisationModel_.ptr());
    surfaceReactionModel_.reset(c.surfaceReactionModel_.ptr());

    dMassDevolatilisation_ = c.dMassDevolatilisation_;
    dMassSurfaceReaction_ = c.dMassSurfaceReaction_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DaemReactingMultiphaseCloud<CloudType>::DaemReactingMultiphaseCloud
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
    daemReactingMultiphaseCloud(),
    cloudCopyPtr_(NULL),
    constProps_(this->particleProperties()),
    devolatilisationModel_(NULL),
    surfaceReactionModel_(NULL),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::DaemReactingMultiphaseCloud<CloudType>::DaemReactingMultiphaseCloud
(
    DaemReactingMultiphaseCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    daemReactingMultiphaseCloud(),
    cloudCopyPtr_(NULL),
    constProps_(c.constProps_),
    devolatilisationModel_(c.devolatilisationModel_->clone()),
    surfaceReactionModel_(c.surfaceReactionModel_->clone()),
    dMassDevolatilisation_(c.dMassDevolatilisation_),
    dMassSurfaceReaction_(c.dMassSurfaceReaction_)
{}


template<class CloudType>
Foam::DaemReactingMultiphaseCloud<CloudType>::DaemReactingMultiphaseCloud
(
    const fvMesh& mesh,
    const word& name,
    const DaemReactingMultiphaseCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    daemReactingMultiphaseCloud(),
    cloudCopyPtr_(NULL),
    constProps_(),
    devolatilisationModel_(NULL),
    surfaceReactionModel_(NULL),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DaemReactingMultiphaseCloud<CloudType>::~DaemReactingMultiphaseCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    label idGas = this->composition().idGas();
    label idLiquid = this->composition().idLiquid();
    label idSolid = this->composition().idSolid();

    parcel.YGas() = this->composition().Y0(idGas);
    parcel.YLiquid() = this->composition().Y0(idLiquid);
    parcel.YSolid() = this->composition().Y0(idSolid);
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    if (fullyDescribed)
    {
        label idGas = this->composition().idGas();
        label idLiquid = this->composition().idLiquid();
        label idSolid = this->composition().idSolid();

        this->checkSuppliedComposition
        (
            parcel.YGas(),
            this->composition().Y0(idGas),
            "YGas"
        );
        this->checkSuppliedComposition
        (
            parcel.YLiquid(),
            this->composition().Y0(idLiquid),
            "YLiquid"
        );
        this->checkSuppliedComposition
        (
            parcel.YSolid(),
            this->composition().Y0(idSolid),
            "YSolid"
        );
    }
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<DaemReactingMultiphaseCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<DaemReactingMultiphaseCloud<CloudType>> td(*this);

        this->solve(td);
    }
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    typedef typename particle::TrackingData<DaemReactingMultiphaseCloud<CloudType>>
        tdType;

    tdType td(*this);

    Cloud<parcelType>::template autoMap<tdType>(td, mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::info()
{
    CloudType::info();

    this->devolatilisation().info(Info);
    this->surfaceReaction().info(Info);
}


template<class CloudType>
void Foam::DaemReactingMultiphaseCloud<CloudType>::writeFields() const
{
    if (this->size())
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


// ************************************************************************* //
