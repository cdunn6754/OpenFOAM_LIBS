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

#include "DaemDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DaemDevolatilisation<CloudType>::
DaemDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    DevolatilisationModel<CloudType>(dict, owner, typeName),
    volatileData_(this->coeffDict().lookup("volatileData")),
    volatileToGasMap_(volatileData_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff")))
{
    if (volatileData_.empty())
    {
        WarningInFunction
            << "Devolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating volatile species:" << endl;

        // Determine mapping between active volatiles and cloud gas components
        const label idGas = owner.composition().idGas(); //thermo index
	// YGasTot: is the initial mass fraction of gas within the particle
        const scalar YGasTot = owner.composition().YMixture0()[idGas];

	// The mass fraction species-wise within the particle
	// this sums to Y_gas not to one
	initialVolatileMassFractions_ = owner.composition().Y0(idGas) * YGasTot;
    }
}


template<class CloudType>
Foam::DaemDevolatilisation<CloudType>::
DaemDevolatilisation
(
    const DaemDevolatilisation<CloudType>& dm
)
:
    DevolatilisationModel<CloudType>(dm),
    volatileData_(dm.volatileData_),
    //YVolatile0_(dm.YVolatile0_),
    volatileToGasMap_(dm.volatileToGasMap_),
    residualCoeff_(dm.residualCoeff_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DaemDevolatilisation<CloudType>::
~DaemDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DaemDevolatilisation<CloudType>::calculate
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
    scalarField& integratedDaemRates
) const
{
    bool done = true;

    // The total current volatile mass fraction of the particle (all species)
    scalar volatileMassFraction = sum(YGasEff);
   
    // If we are restarting, there is no disk storage for "done" so it will try to 
    // devolatilise parcels which are already done. This leads to dividing by zero 
    // in the code below since volatileMassFraction == 0. in this case. So probably 
    // we should not let that happen
    if (volatileMassFraction >= VSMALL)
      {

	// Quadrature points and wieghts
	scalarField quadPoints(4,0);
	scalarField quadWeights(4,0);
	quadPoints = volatileData_[0].x();
	quadWeights = volatileData_[0].W();

        // Model coefficients
        const scalar A0 = volatileData_[0].A0();
        const scalar Emean = volatileData_[0].Emean();
	const scalar Estd = volatileData_[0].Estd();
	const scalar m = volatileData_[0].m();
	
	// Determine the sample Activation Energies from the quad points
	scalarField Esamp = Emean + quadPoints * pow(2., 0.5) * Estd * m;

	// V_star is the initial mass fraction of all volatile species in the coal
	// or alternatively V* is V as t -> infty where V is the massfrac evolved at t
	scalar V_star = sum(initialVolatileMassFractions_);

        // taking a sum over all quad points for this timestep, should start at 0
	scalar sum = 0.;

	// loop through all of the quadrature points 
	forAll(integratedDaemRates, j)
	  {
	    // Define the quantities for this iteration (this quad point)
	    scalar Ei = Esamp[j];
	    scalar Wi = quadWeights[j];
	    
            // Integrate by this particle time step dt
	    integratedDaemRates[j] = integratedDaemRates[j]
	      + dt * A0 * exp(-Ei/(T*(RR/1000.0)));
	    
	    // Sum over all of the quadrature points
	    sum = sum + Wi * (m/pow(3.1415, 0.5)) * 
	      pow(2.718281, -pow(Ei - Emean, 2.)/(2.* pow(Estd, 2.)))
	      * pow(2.718281, -integratedDaemRates[j]);
	  }

	scalar V = V_star * (1. - sum); // update the current volatile mass fraction

	scalar particleNewMass = (1.0 - V)*mass0;

	Info << "New Mass " << mass0 - particleNewMass << endl;

	// Every species will have its own volatileData object,
	// volatileData_ stores them, so iterate through them here
	// and update the mass loss on a species-wise basis
	forAll(initialVolatileMassFractions_, i)
	  {
	    //masses for this species
	    const scalar massVolatile0 = mass0*initialVolatileMassFractions_[i];
	    const scalar massVolatile = mass*YGasEff[i];

	    
	    // The mass fraction of this species within the gas mass
	    scalar speciesMassFraction = YGasEff[i]/volatileMassFraction;

	    // The mass of this species lost during this iteration
	    scalar speciesMassLoss = speciesMassFraction * (mass - particleNewMass);

	    // During manual restart (e.g. setting all to IDR = 0.0)
	    // can have problems with new particle mass being more than current mass
	    // and predicting negative mass loss. This happens when they are 
	    // nearly out of volatile matter anyways, so use the constant rate
	    // method from the OF devolatilisation model with A0=12
	    if (speciesMassLoss <= 0.0)
	      {
		// Mass transfered from the particle to the carrier gas phase
		// on a per species basis
		dMassDV[i] = min(massVolatile, dt * 12. * massVolatile0);
		Info << "WARNING: Forced to use constant rate model " << endl;
	      }
	    else
	      {
		// Mass transfered from the particle to the carrier gas phase
		// on a per species basis
		dMassDV[i] = min(massVolatile, speciesMassLoss);
	      }


	    
	    // Combustion allowed once all volatile components evolved
	    done = done && (massVolatile <= residualCoeff_*massVolatile0);
	  }
      }

    else // if there are not volatiles to begin with
      {
	dMassDV[0] = 0.0;
	dMassDV[1] = 0.0;
	dMassDV[2] = 0.0;
      }

    Info << sum(dMassDV)  << endl;

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
