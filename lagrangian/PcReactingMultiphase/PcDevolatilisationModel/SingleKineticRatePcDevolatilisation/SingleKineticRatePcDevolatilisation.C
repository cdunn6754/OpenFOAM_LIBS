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

	// Setting the proper value for Ydaf0
	const label idSolid = owner.composition().idSolid();
	const label idLiquid = owner.composition().idLiquid();
	const scalarField& YSolid0 = owner.composition().Y0(idSolid);
	const scalarField& YLiquid0 = owner.composition().Y0(idLiquid);
	const label ashId = owner.composition().localId(idSolid, "ash");
	const label waterId = owner.composition().localId(idLiquid, "H2O");
	const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
	const scalar YLiquidTot = owner.composition().YMixture0()[idLiquid];
	const scalar Yash = YSolidTot * YSolid0[ashId];
	const scalar Ywater = YLiquidTot * YLiquid0[waterId];
	Ydaf0_ = 1.0 - Yash - Ywater;
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
    Info << "\nTar Fields: " << tarFields[0] << endl;

    // Initial daf mass of particle
    const scalar dafMass0 = this->Ydaf0_ * mass0;

    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];
        const scalar massVolatile0 = mass0*YVolatile0_[i];
        const scalar massVolatile = mass*YGasEff[id];

	// For the PC coal lab devol rate laws we need daf based 
	// YdafVolatile0, as opposed to the YVolatile0 we already have.
	const scalar YdafVolatile0 = YVolatile0_[i]/this->Ydaf0_;

	// Find the percentage of the intial daf mass 
	// that has been devolatilized
	// TDevoled starts at 0.0 and approaches YdafVolatile0
	const scalar massDevoled = massVolatile0 - massVolatile;
	scalar YDevoled = massDevoled/dafMass0;	
 

        // Combustion allowed once all volatile components evolved
        done = done && (massVolatile <= residualCoeff_*massVolatile0);

        // Model coefficients
        const scalar Ap = volatileData_[i].Ap();
        const scalar Ep = volatileData_[i].Ep();

	// make sure we dont exceed the amount of volatiles in the coal
	if (YDevoled >= YdafVolatile0)
	  {
	    YDevoled = YdafVolatile0;
	  }

	
	// Convert the builtin RR from [J/kmol K] 
	// to [kcal/mol K]
	const scalar pcR = ((RR/1000.)/4184);

        // Kinetic rate
        const scalar kappa = Ap*exp(-Ep/(pcR*T));

        // Mass transferred from particle to carrier gas phase
	// The rates are also based on the daf mass
	const scalar massTransfered = dt * kappa * (YdafVolatile0 - YDevoled) 
	  * dafMass0;
        dMassDV[id] = min(massTransfered, massVolatile);
    }

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
