/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "SingleKineticRatePcDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SingleKineticRatePcDevolatilisation<CloudType>::
SingleKineticRatePcDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PcDevolatilisationModel<CloudType>(dict, owner, typeName),
    volatileData_(this->coeffDict().lookup("volatileData")),
    YVolatile0_(volatileData_.size()),
    volatileToGasMap_(volatileData_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff"))),
    Ydaf0_(1.0)
{
    if (volatileData_.empty())
    {
        WarningInFunction
            << "PcDevolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating volatile species:" << endl;

        // Determine mapping between active volatiles and cloud gas components
        const label idGas = owner.composition().idGas();
        const scalar YGasTot = owner.composition().YMixture0()[idGas];
        const scalarField& YGas0 = owner.composition().Y0(idGas);
        forAll(volatileData_, i)
        {
            const word& specieName = volatileData_[i].name();
            const label id = owner.composition().localId(idGas, specieName);
            volatileToGasMap_[i] = id;
            YVolatile0_[i] = YGasTot*YGas0[id];

            Info<< "    " << specieName << ": particle mass fraction = "
                << YVolatile0_[i] << endl;
        }

	// Setting the proper value for Ydaf0_
	const label ashId = owner.composition().localId(idGas, "ash");
	const label waterId = owner.composition().localId(idGas, "H2O");
	Ydaf0_ = 1.0 - YGas0[ashId] - YGas0[waterId];
    }
}


template<class CloudType>
Foam::SingleKineticRatePcDevolatilisation<CloudType>::
SingleKineticRatePcDevolatilisation
(
    const SingleKineticRatePcDevolatilisation<CloudType>& dm
)
:
    PcDevolatilisationModel<CloudType>(dm),
    volatileData_(dm.volatileData_),
    YVolatile0_(dm.YVolatile0_),
    volatileToGasMap_(dm.volatileToGasMap_),
    residualCoeff_(dm.residualCoeff_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SingleKineticRatePcDevolatilisation<CloudType>::
~SingleKineticRatePcDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SingleKineticRatePcDevolatilisation<CloudType>::calculate
(
    const scalar dt,
    const scalar age,
    const scalar mass0,
    const scalar mass,
    const scalar T,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalarField& tarFields
) const
{
    bool done = true;

    // Initial daf mass of particle
    const scalar dafMass0 = this->Ydaf0_ * mass0;

    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];
        const scalar massVolatile0 = mass0*YVolatile0_[i];
        const scalar massVolatile = mass*YGasEff[id];

	const scalar massDevoled = massVolatile0 - massVolatile;
	scalar fracDevoled = massDevoled/dafMass0;
	
 

        // Combustion allowed once all volatile components evolved
        done = done && (massVolatile <= residualCoeff_*massVolatile0);

        // Model coefficients
        const scalar A1 = volatileData_[i].A1();
        const scalar E = volatileData_[i].E();
	const scalar Y_inf = volatileData_[i].YVolinf();

	// make sure we dont exceed the amount of volatiles in the coal
	if (fracDevoled >= Y_inf)
	  {
	    fracDevoled = Y_inf;
	  }

        // Kinetic rate
        const scalar kappa = A1*exp(-E/(RR*T));

        // Mass transferred from particle to carrier gas phase
	const scalar massTransfered = dt*kappa*(Y_inf - fracDevoled) * dafMass0;
        dMassDV[id] = min(massTransfered, massVolatile);
    }

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
