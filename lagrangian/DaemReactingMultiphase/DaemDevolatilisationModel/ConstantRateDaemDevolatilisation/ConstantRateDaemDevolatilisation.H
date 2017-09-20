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
    Foam::ConstantRateDaemDevolatilisation

Description
    Constant rate devolatisation model
    - need to set vapourisation temperature to 600 K

\*---------------------------------------------------------------------------*/

#ifndef ConstantRateDaemDevolatilisation_H
#define ConstantRateDaemDevolatilisation_H

#include "DaemDevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                 Class ConstantRateDaemDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConstantRateDaemDevolatilisation
:
    public DaemDevolatilisationModel<CloudType>
{
    // Private data

        // Model constants

            //- List of volatile data - (name A0)
            List<Tuple2<word, scalar>> volatileData_;

            //- List of initial volatile mass fractions
            List<scalar> YVolatile0_;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;


public:

    //- Runtime type information
    TypeName("constantRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        ConstantRateDaemDevolatilisation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        ConstantRateDaemDevolatilisation
        (
            const ConstantRateDaemDevolatilisation<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DaemDevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<DaemDevolatilisationModel<CloudType>>
            (
                new ConstantRateDaemDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ConstantRateDaemDevolatilisation();


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
    #include "ConstantRateDaemDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
