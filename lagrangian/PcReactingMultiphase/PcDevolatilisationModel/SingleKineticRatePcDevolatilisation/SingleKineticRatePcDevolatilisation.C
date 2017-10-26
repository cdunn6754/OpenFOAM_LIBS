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
    scalarField& tarProps,
    scalarField& dMassSP
) const
{
    bool done = true;

    // Initial daf mass of particle
    // needs to incorporate the tar mass 
    // that isnt included
    // so we are basically tacking on a fake 
    // mass here
    const scalar dafMass0 = Ydaf0_ * mass0;

    // Convert the builtin RR from [J/kmol K] 
    // to [kcal/mol K] (PCCL B.S.)
    const scalar pcR = ((RR/1000.)/4184);

    // Find the current daf mass fraction of tar within the particle
    //const scalar Yash = YSolidEff[ashId_];
    //const scalar Ywater = YLiquidEff[waterId_];
    //const scalar YdafTar = (YdafInfTar_ - tarFields[0]);
    //const scalar dafMass = (1.0 - Ywater - Yash) * mass;

    


    // Loop through volatile species
    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];

	// Find the mass of this specie within the particle
        const scalar massVolatile0 = mass0*YVolatile0_[i];
	// get the current mass of this specie volatiles in the particle
        const scalar massVolatile = mass*YGasEff[id];

	// For the PC coal lab devol rate laws we need daf based
	// YdafVolatile0, as opposed to the YVolatile0 we already have.
	const scalar YdafVolatile0 = (YVolatile0_[i]/Ydaf0_);

	// YdafDevoled is the mass fraction of current specie mass lost
	// compared to the initial daf mass (dafMass0)
	// YdafDevoled starts at 0.0 and approaches YdafVolatile0
	const scalar massDevoled = massVolatile0 - massVolatile;
	scalar YdafDevoled = massDevoled/(dafMass0);	
 
	Info << "YdafDevoled " <<  volatileData_[i].name() << ": "
	     << YdafDevoled << endl;

	Info << "YdafDevol0 " <<  volatileData_[i].name() << ": "
	     << YdafVolatile0 << endl;

        // Combustion allowed once all volatile components evolved
        done = done && (massVolatile <= residualCoeff_*massVolatile0);

        // Model coefficients
        const scalar Ap = volatileData_[i].Ap();
        const scalar Ep = volatileData_[i].Ep();

	// make sure we dont exceed the amount of volatiles in the coal
	if (YdafDevoled >= YdafVolatile0)
	  {
	    YdafDevoled = YdafVolatile0;
	  }

        // Kinetic rate
        const scalar kappa = Ap*exp(-Ep/(pcR*T));

        // Mass transferred from particle to carrier gas phase
	// The rates are also based on the daf mass
	const scalar YdafTransfered = dt * kappa * (YdafVolatile0 - YdafDevoled);


	// Give the primary devolatilization products to gas phase
	// (including tar, we will remove in outer function)
	dMassDV[id] = min(YdafTransfered*dafMass0, massVolatile);

	// In the case of Tar we need to do the secondary pyrolysis
	// in addition to the primary devolatilization
	if (volatileData_[i].name() == "TAR")
	  {
	    // - Get the particle data from tarProps
	    // Primary released by particle daf mass fraction
	    scalar& Yp = tarProps[0];
	    // difference between Primary tar released (Yp)
	    // and actual tar remaining (i.e. secondary tar (Ys))
	    scalar& dY = tarProps[1];

	    // increment the primary tar release
	    Yp += YdafTransfered;

	    // - Handle the Secondary pyrolysis
	    // This corresponds to EQN 9 chapter 13 of PCCL book
	    // The primary tar is broken down

	    const scalar As = volatileData_[i].As();
	    const scalar Es = volatileData_[i].Es();
	    // get the daf mass fraction of tar within the particle initially
	    const scalar Yinf = volatileData_[i].Yinfs(); 

	    const scalar kappa = As * exp(-Es/(pcR*T));
	    const scalar dyIncrement = dt * kappa * (Yinf - dY);

	    // increment the dY value of this particle
	    dY += dyIncrement;
	    
	    // mass of this mass fraction increment
	    const scalar massIncrement = dyIncrement * dafMass0;

	    // Separate the tar breakdown mass into its contituents
	    // and release them to the gas phase
	    forAll(volatileData_, j)
	      {
		if (volatileData_[j].name() == "TAR")
		  {
		    continue;
		  }
		const label ids = volatileToGasMap_[j];
		dMassSP[ids] = volatileData_[j].Yinfs() * (massIncrement);
	      }
	    
	    Info << "Yp: " << Yp << endl;
	    Info << "dY: " << dY << endl;
	  } // end tar

    } // end per specie loop



    // // --- Secondary Pyrolysis
    // // - Get the particle data from tarPro
    // // Primary released by particle daf mass fraction
    // scalar& Yp = tarFields[0];
    // // difference between Primary tar released (Yp)
    // // and actual tar remaining (i.e. secondary tar (Ys))
    // scalar& dY = tarFields[1];


    // // - Devolatilize the primary tar (just like any other species)
    // const scalar kappaTar = Atar_ * exp(- Etar_/(pcR*T));
    // const scalar YpNew = dt * kappaTar * (YdafInfTar_ - Yp);
     
    // // Add to primary tar and decompose tar mass fraction
    // Yp += YpNew;


    // // -  Now evaluate the dY equation (Y here based on daf)

    // // rate constant for dY equation
    // const scalar kappaDy = AtarDy_ * exp(- EtarDy_/(pcR*T) );
    // const scalar dYnew = dt * kappaDy * (YdafInfTar_ - dY);
    // dY += dYnew;
    // const scalar dYmass = dYnew * dafMass0;

    // // Separate the tar breakdown mass into its contituents
    // forAll(volatileData_, j)
    //   {
    // 	const label ids = volatileToGasMap_[j];
    // 	dMassDV[ids] += volatileData_[j].Yinfs() * (dYmass);
    //   }

    if (done && canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
