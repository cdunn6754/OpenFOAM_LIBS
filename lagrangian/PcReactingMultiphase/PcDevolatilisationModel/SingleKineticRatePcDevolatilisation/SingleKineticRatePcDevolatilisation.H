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

Class
    Foam::SingleKineticRatePcDevolatilisation

Description
    Single kinetic rate devolatisation model for the 
    primary particle pyrolysis (i.e. devolatilization). This
    class also models secondary pyrolysis.

    Primary Pyrolysis:
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = A*exp(- E/R*T)

    Where:
        kappa = rate constant
        A    = activation energy (user input) [1/s]
        E     = pre-exponential factor (user input) [kcal/mole]
        R     = universal gas constant
        T     = temperature


    Secondary Pyrolysis
    - speciated as well
    - presumes that one of the species of primary pyrolysis is
      named "TAR". This tar is not released to the gas phase
      as the other species are and is instead accumulated in the
      parcel. 
    - Three variables related to the secondary pyrolysis are tracked
      with the particle. Yp, dY and Ysoot.
     . Yp is the product of the primary pyrolysis, i.e. it is the amount
       of TAR that has been devolatilized (only continues to increase)
     . dY is the difference between Yp and Ys. Where Ys is the actual
       amount of tar remaining in the particle (i.e. Yp - (the amount
       that has been lost to the gas phase through secondary pyrolysis)). 
       Yp >= Ys and once primary pyrolysis is complete (=> Yp=Yinf) and 
       all of the tar has been pyrolyzed to gas phase species 
       => dY=Yp=Yinf. (where Yinf is the daf mass fraction of tar in the
       parcel.
     . Ysoot  TBD, currently just released as an inert gas species


    Usage:

        SingleKineticRateDevolatilisationCoeffs
        {
            volatileData
            (
                (CH4     12     0.5    191.0    8.3   .23) 
                (CO2     12     0.5    191.0    8.3   .35)
            );

            volatileResidualCoeff 1e-6;
        }

	The parameters are passed in in this form
	(speciesName    Ap    Ep    As    Ap    YInfs)

	Here p refers to primary rates, s refers to secondary 
	pyrolyis (tar breakdown). Each species that is involved in
	either process will need a complete entry. Their inactivity 
	in the other process can be ensured by setting the activation
	energy to zero. See units above near discussion of kappa.

	YInfs is the secondary pyrolysis ultimate yield. This should be
	the mass fraction of the tar which this species comprises. A primary 
	pyrolysis counterpart does not need to be specified here as it is read 
	from the Cloud composition. 
	

\*---------------------------------------------------------------------------*/

#ifndef SingleKineticRatePcDevolatilisation_H
#define SingleKineticRatePcDevolatilisation_H

#include "PcDevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class SingleKineticRatePcDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SingleKineticRatePcDevolatilisation
:
    public PcDevolatilisationModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volatileData
    {
        // Private data

            //- Specie name
            word name_;

            //- Activation energy Primary
            scalar Ap_;

            //- Pre-exponential factor Primary
            scalar Ep_;

            //- Activation energy Secondary
            scalar As_;

            //- Pre-exponential factor Secondary
            scalar Es_;

            //- Ultimate species mass fraction (within tar)
            scalar Yinfs_;

    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null),
                Ap_(0.0),
                Ep_(0.0),
                As_(0.0),
                Es_(0.0),
		Yinfs_(0.0)
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is),
                Ap_(readScalar(is)),
                Ep_(readScalar(is)),
                As_(readScalar(is)),
                Es_(readScalar(is)),
                Yinfs_(readScalar(is))
            {}

            //- Construct as copy
            volatileData(const volatileData& vd)
            :
                name_(vd.name_),
                Ap_(vd.Ap_),
                Ep_(vd.Ep_),
                As_(vd.As_),
                Es_(vd.Es_),
                Yinfs_(vd.Yinfs_)
            {}


        //- Destructor
        ~volatileData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the activation energy
                scalar Ap() const
                {
                    return Ap_;
                }

                //- Return const access to the pre-exponential factor
                scalar Ep() const
                {
                    return Ep_;
                }
                //- Return const access to the activation energy
                scalar As() const
                {
                    return As_;
                }

                //- Return const access to the pre-exponential factor
                scalar Es() const
                {
                    return Es_;
                }
                //- What do you suppose this function does?
                scalar Yinfs() const
                {
                    return Yinfs_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.Ap_  >> vd.Ep_ >> 
		  vd.As_ >> vd.Es_ >> vd.Yinfs_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.Ap_ << token::SPACE
                    << vd.Ep_ << token::SPACE
                    << vd.As_ << token::SPACE
                    << vd.Es_ << token::SPACE
		    << vd.Yinfs_
                    << token::END_LIST;

                return os;
            }
    };


    // Private data

        // Model constants

            //- List of volatile data - (name A1 E)
            List<volatileData> volatileData_;

            //- List of initial volatile mass fractions
            List<scalar> YVolatile0_;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;

  // In order to work with daf basis in main funtion
  // we need the mass fraction of ash and water of
  // the particles so Ydaf0 is the initial 
  // mass fraction of the particle that is daf
  // it is assumed that these fractions are the same for all particles
  scalar Ydaf0_;

  // ashID and waterID for the Composition object
  label ashId_;
  label waterId_;

  //- Ultimate daf yield of tar
  scalar YdafInfTar_;
  //- Activatation Energy for primary tar equation [1/s]
  scalar Atar_;
  //- Frequency Factor for primary tar equation [kcal/mol]
  scalar Etar_;
  //- Activatation Energy for dY tar equation [1/s]
  scalar AtarDy_;
  //- Frequency Factor for dY tar equation [kcal/mol]
  scalar EtarDy_;

  //- Reference to the cloud
  const CloudType& owner_;


public:

    //- Runtime type information
    TypeName("singleKineticRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        SingleKineticRatePcDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        SingleKineticRatePcDevolatilisation
        (
            const SingleKineticRatePcDevolatilisation<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<PcDevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<PcDevolatilisationModel<CloudType>>
            (
                new SingleKineticRatePcDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SingleKineticRatePcDevolatilisation();


    // Member Functions

        //- Update model
        virtual void calculate
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
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SingleKineticRatePcDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
