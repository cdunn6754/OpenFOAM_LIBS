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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::NotionalParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties()
{}


template<class ParcelType>
inline Foam::NotionalParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp)
{}


template<class ParcelType>
inline Foam::NotionalParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict)
{}


template<class ParcelType>
inline Foam::NotionalParcel<ParcelType>::NotionalParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPtI
)
:
    ParcelType(mesh, position, celli, tetFacei, tetPtI)
    // YGas_(0),
    // YLiquid_(0),
    // YSolid_(0),
    // canCombust_(0)
{}


template<class ParcelType>
inline Foam::NotionalParcel<ParcelType>::NotionalParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalarField& Y0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        position,
        celli,
        tetFacei,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
	YGas0,
	YLiquid0,
	YSolid0,
        constProps
    )
    // YGas_(YGas0),
    // YLiquid_(YLiquid0),
    // YSolid_(YSolid0),
    // canCombust_(0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

// clint 07-17-17 not needed anymore they are already defined in RMparcel

// template<class ParcelType>
// inline Foam::scalar
// Foam::NotionalParcel<ParcelType>::constantProperties::TDevol() const
// {
//     return TDevol_.value();
// }


// template<class ParcelType>
// inline Foam::scalar
// Foam::NotionalParcel<ParcelType>::constantProperties::LDevol() const
// {
//     return LDevol_.value();
// }


// template<class ParcelType>
// inline Foam::scalar
// Foam::NotionalParcel<ParcelType>::constantProperties::
// hRetentionCoeff() const
// {
//     scalar value = hRetentionCoeff_.value();

//     if ((value < 0) || (value > 1))
//     {
//         FatalErrorInFunction
//             << "hRetentionCoeff must be in the range 0 to 1" << nl
//             << exit(FatalError) << endl;
//     }

//     return value;
// }


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

// template<class ParcelType>
// inline const Foam::scalarField& Foam::NotionalParcel<ParcelType>::
// YGas() const
// {
//     return YGas_;
// }


// template<class ParcelType>
// inline const Foam::scalarField& Foam::NotionalParcel<ParcelType>::
// YLiquid() const
// {
//     return YLiquid_;
// }


// template<class ParcelType>
// inline const Foam::scalarField& Foam::NotionalParcel<ParcelType>::
// YSolid() const
// {
//     return YSolid_;
// }


// template<class ParcelType>
// inline Foam::label
// Foam::NotionalParcel<ParcelType>::canCombust() const
// {
//     return canCombust_;
// }


// template<class ParcelType>
// inline Foam::scalarField& Foam::NotionalParcel<ParcelType>::YGas()
// {
//     return YGas_;
// }


// template<class ParcelType>
// inline Foam::scalarField& Foam::NotionalParcel<ParcelType>::YLiquid()
// {
//     return YLiquid_;
// }


// template<class ParcelType>
// inline Foam::scalarField& Foam::NotionalParcel<ParcelType>::YSolid()
// {
//     return YSolid_;
// }


// template<class ParcelType>
// inline Foam::label& Foam::NotionalParcel<ParcelType>::canCombust()
// {
//     return canCombust_;
// }


// ************************************************************************* //
