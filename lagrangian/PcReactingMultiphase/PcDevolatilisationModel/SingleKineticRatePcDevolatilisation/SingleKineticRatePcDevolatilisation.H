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
    Single kinetic rate devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = A1.exp(- E/R.T)

    Where:
        kappa = rate constant
        A1    = activation energy (user input)
        E     = pre-exponential factor (user input)
	Yinf  = Total mass fraction of tar in daf basis
        R     = universal gas constant
        T     = temperature

    Usage:

        SingleKineticRatePcDevolatilisationCoeffs
        {
            volatileData
            (
                (CH4     12     0.5    0.03)   // (name A1 E Yinf)
                (CO2     12     0.5    0.03)   // (name A1 E Yinf)
            );

            volatileResidualCoeff 1e-6;
        }

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

            //- Activation energy
            scalar A1_;

            //- Pre-exponential factor
            scalar E_;

      // Ultimate Yield in daf basis 
      // (ie what fraction of daf initial mass is volatiles? 
      // YVolinf_, thats what. 
      // Then it is also the total fraction that can be 
      // devolatilized eventually (Clint 10-17)
      scalar YVolinf_;

    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null),
                A1_(0.0),
                E_(0.0),
		YVolinf_(0.0)
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is),
                A1_(readScalar(is)),
                E_(readScalar(is)),
                YVolinf_(readScalar(is))
            {}

            //- Construct as copy
            volatileData(const volatileData& vd)
            :
                name_(vd.name_),
                A1_(vd.A1_),
                E_(vd.E_),
                YVolinf_(vd.YVolinf_)
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
                scalar A1() const
                {
                    return A1_;
                }

                //- Return const access to the pre-exponential factor
                scalar E() const
                {
                    return E_;
                }

                scalar YVolinf() const
                {
                    return YVolinf_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.A1_  >> vd.E_ >> vd.YVolinf_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.A1_ << token::SPACE
                    << vd.E_ << token::SPACE
		    << vd.YVolinf_
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
