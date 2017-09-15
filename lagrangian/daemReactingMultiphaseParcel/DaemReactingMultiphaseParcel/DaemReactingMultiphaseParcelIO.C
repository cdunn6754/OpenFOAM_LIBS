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

#include "DaemReactingMultiphaseParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::DaemReactingMultiphaseParcel<ParcelType>::propertyList_ =
    Foam::DaemReactingMultiphaseParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::DaemReactingMultiphaseParcel<ParcelType>::sizeofFields_
(
    0
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DaemReactingMultiphaseParcel<ParcelType>::DaemReactingMultiphaseParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{
  //Info << "\nConstructed from ISTREAM \n" << endl;
  
    if (readFields)
    {
        DynamicList<scalar> Yg;
        DynamicList<scalar> Yl;
        DynamicList<scalar> Ys;
	DynamicList<scalar> IDR;


        is >> Yg >> Yl >> Ys >> IDR;

        YGas_.transfer(Yg);
        YLiquid_.transfer(Yl);
        YSolid_.transfer(Ys);
	integratedDaemRates_.transfer(IDR);

        // scale the mass fractions
        const scalarField& YMix = this->Y_;
        YGas_ /= YMix[GAS] + ROOTVSMALL;
        YLiquid_ /= YMix[LIQ] + ROOTVSMALL;
        YSolid_ /= YMix[SLD] + ROOTVSMALL;
    }

    // Check state of Istream
    is.check
    (
        "DaemReactingMultiphaseParcel<ParcelType>::DaemReactingMultiphaseParcel"
        "("
            "const polyMesh&, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::DaemReactingMultiphaseParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::DaemReactingMultiphaseParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
  Info << "\nIn read fields" << endl;
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c, compModel);

    // Get names and sizes for each Y...
    const label idGas = compModel.idGas();
    const wordList& gasNames = compModel.componentNames(idGas);
    const label idLiquid = compModel.idLiquid();
    const wordList& liquidNames = compModel.componentNames(idLiquid);
    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);
    const wordList& stateLabels = compModel.stateLabels();

    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<DaemReactingMultiphaseParcel<ParcelType>>, c, iter)
    {
        DaemReactingMultiphaseParcel<ParcelType>& p = iter();
        p.YGas_.setSize(gasNames.size(), 0.0);
        p.YLiquid_.setSize(liquidNames.size(), 0.0);
        p.YSolid_.setSize(solidNames.size(), 0.0);
    }

    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        IOField<scalar> YGas
        (
            c.fieldIOobject
            (
                "Y" + gasNames[j] + stateLabels[idGas],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<DaemReactingMultiphaseParcel<ParcelType>>,
            c,
            iter
        )
        {
            DaemReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YGas_[j] = YGas[i++]/(p.Y()[GAS] + ROOTVSMALL);
        }
    }
    // Populate YLiquid for each parcel
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid
        (
            c.fieldIOobject
            (
                "Y" + liquidNames[j] + stateLabels[idLiquid],
                 IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<DaemReactingMultiphaseParcel<ParcelType>>,
            c,
            iter
        )
        {
            DaemReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YLiquid_[j] = YLiquid[i++]/(p.Y()[LIQ] + ROOTVSMALL);
        }
    }
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j] + stateLabels[idSolid],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<DaemReactingMultiphaseParcel<ParcelType>>,
            c,
            iter
        )
        {
            DaemReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YSolid_[j] = YSolid[i++]/(p.Y()[SLD] + ROOTVSMALL);
        }
    }

    // DAEM read stuff 06-17 Clint
    IOField<scalarField> 
      integratedDaemRates(c.fieldIOobject("integratedDaemRates",
					  IOobject::MUST_READ));
    c.checkFieldIOobject(c, integratedDaemRates);

    label i = 0;
    forAllIter(typename CloudType, c, iter)
      {
	DaemReactingMultiphaseParcel<ParcelType>& p = iter();
	
	p.integratedDaemRates_ = integratedDaemRates[i];

	i++;
      }
}


template<class ParcelType>
template<class CloudType>
void Foam::DaemReactingMultiphaseParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    //label np=c.size();
    

}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::DaemReactingMultiphaseParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    label np = c.size();

    // Write the composition fractions
    if (np > 0)
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            IOField<scalar> YGas
            (
                c.fieldIOobject
                (
                    "Y" + gasNames[j] + stateLabels[idGas],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<DaemReactingMultiphaseParcel<ParcelType>>,
                c,
                iter
            )
            {
                const DaemReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YGas[i++] = p0.YGas()[j]*p0.Y()[GAS];
            }

            YGas.write();
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid
            (
                c.fieldIOobject
                (
                    "Y" + liquidNames[j] + stateLabels[idLiquid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<DaemReactingMultiphaseParcel<ParcelType>>,
                c,
                iter
            )
            {
                const DaemReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YLiquid[i++] = p0.YLiquid()[j]*p0.Y()[LIQ];
            }

            YLiquid.write();
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid
            (
                c.fieldIOobject
                (
                    "Y" + solidNames[j] + stateLabels[idSolid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<DaemReactingMultiphaseParcel<ParcelType>>,
                c,
                iter
            )
            {
                const DaemReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YSolid[i++] = p0.YSolid()[j]*p0.Y()[SLD];
            }

            YSolid.write();
        }
    }


    // DAEM rates writeout added 06-17 Clint 
    IOField<scalarField> 
      integratedDaemRates(c.fieldIOobject("integratedDaemRates", 
					  IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
      {
	const DaemReactingMultiphaseParcel<ParcelType>& p = iter();
	integratedDaemRates[i] = p.integratedDaemRates();

	i++;
      }

    integratedDaemRates.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DaemReactingMultiphaseParcel<ParcelType>& p
)
{
    scalarField YGasLoc(p.YGas()*p.Y()[0]);
    scalarField YLiquidLoc(p.YLiquid()*p.Y()[1]);
    scalarField YSolidLoc(p.YSolid()*p.Y()[2]);
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << YGasLoc
            << token::SPACE << YLiquidLoc
            << token::SPACE << YSolidLoc
	    << token::SPACE << p.integratedDaemRates();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os  << YGasLoc << YLiquidLoc << YSolidLoc 
	    << p.integratedDaemRates();
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const DaemReactingMultiphaseParcel<ParcelType>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
