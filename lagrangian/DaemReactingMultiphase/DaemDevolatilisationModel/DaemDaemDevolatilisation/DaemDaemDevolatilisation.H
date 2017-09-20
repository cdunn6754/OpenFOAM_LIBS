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
    Foam::DaemDaemDevolatilisation

Description
    DAEM  devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = A1.exp(- E/R.T)

    Where:
        kappa = rate constant
        A1    = activation energy (user input)
        E     = pre-exponential factor (user input)
        R     = universal gas constant
        T     = temperature

    Usage:

        DaemDaemDevolatilisationCoeffs
        {
            volatileData
            (
                (CH4     12     0.5)   // (name A1 E)
                (CO2     12     0.5)   // (name A1 E)
            );

            volatileResidualCoeff 1e-6;
        }

\*---------------------------------------------------------------------------*/

#ifndef DaemDaemDevolatilisation_H
#define DaemDaemDevolatilisation_H

#include "DaemDevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class DaemDaemDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DaemDaemDevolatilisation
:
    public DaemDevolatilisationModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volatileData
    {
      // Private data

      //- Specie name
      word name_;

      //- Pre-exponential
      scalar A0_;

      //- Activation Energy mean
      scalar Emean_;

      //- Activation Energy standard deviation
      scalar Estd_;

      // Gaussian-Hermite Quadrature points has four scalar values
      scalarField x_;
      
      // Gaussian-Hermite Quadrature weights has four scalar values
      scalarField W_;
      
      // Quadrature Scaling factor
      scalar m_;
      


    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
	      name_(word::null)
	      //W_(scalarField temp(4,0))
                // A1_(0.0),
                // E1_(0.0)
            {
	      scalarField temp(4,0.);
	      W_ = temp;
	      x_ = temp;
	    }

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is)
                // A1_(readScalar(is)),
                // E1_(readScalar(is)),
	    	// E2_(readScalar(is))
            {
	    }

            //- Construct as copy
            volatileData(const volatileData& vd)
            :
                name_(vd.name_),
                A0_(vd.A0_),
                Emean_(vd.Emean_),
		Estd_(vd.Estd_),
		x_(vd.x_),
		W_(vd.W_),
		m_(vd.m_)

            {
	    }


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

                //- Return const access to the pre-exponential factor
                scalar A0() const
                {
                    return A0_;
                }

                //- Return const access to the activation energy mean
                scalar Emean() const
                {
                    return Emean_;
                }

                //- Return const access to the activation energy std
                scalar Estd() const
                {
                    return Estd_;
                }

                //- Return const access to the quad points
                scalarField x() const
                {
                    return x_;
                }

                //- Return const access to the third quad weight
                scalarField W() const
                {
		  return W_;
		}

                //- Return const access to the quad scaling factor
                scalar m() const
                {
                    return m_;
                }



        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.A0_  >> vd.Emean_ >> vd.Estd_ >>
		  vd.x_ >>
		  vd.W_ >> vd.m_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.A1_ << token::SPACE
                    << vd.E1_
                    << token::END_LIST;

                return os;
            }
    };


    // Private data

        // Model constants

            //- List of volatile data objects, one for each species
            List<volatileData> volatileData_;

  //- List of initial volatile mass fractions for each species
  // does not sum to one. It is the mass fraction wrt the total particle mass
  // not the gas fraction of the particle
  scalarField initialVolatileMassFractions_;


            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;

public:

    //- Runtime type information
    TypeName("DaemDevolatilisation");


    // Constructors

        //- Construct from dictionary
        DaemDaemDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        DaemDaemDevolatilisation
        (
            const DaemDaemDevolatilisation<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DaemDevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<DaemDevolatilisationModel<CloudType>>
            (
                new DaemDaemDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DaemDaemDevolatilisation();


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
	    scalarField& integratedDaemRates
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DaemDaemDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
