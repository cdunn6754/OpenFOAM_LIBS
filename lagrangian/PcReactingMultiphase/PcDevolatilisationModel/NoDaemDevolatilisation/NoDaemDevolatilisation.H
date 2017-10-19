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
    Foam::NoDaemDevolatilisation

Description
    Dummy devolatilisation model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoDaemDevolatilisation_H
#define NoDaemDevolatilisation_H

#include "DaemDevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class NoDaemDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoDaemDevolatilisation
:
    public DaemDevolatilisationModel<CloudType>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoDaemDevolatilisation(const dictionary&, CloudType& owner);

        //- Construct copy
        NoDaemDevolatilisation(const NoDaemDevolatilisation<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DaemDevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<DaemDevolatilisationModel<CloudType>>
            (
                new NoDaemDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoDaemDevolatilisation();


    // Member Functions

        //- Flag to indicate whether model activates devolatilisation model
        virtual bool active() const;

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
    #include "NoDaemDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
