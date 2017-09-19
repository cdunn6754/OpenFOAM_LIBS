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
    Foam::DaemDevolatilisationModel

Description
    Templated devolatilisation model class

SourceFiles
    DaemDevolatilisationModel.C
    DaemDevolatilisationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef DaemDevolatilisationModel_H
#define DaemDevolatilisationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DaemDevolatilisationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DaemDevolatilisationModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //- Runtime type information
    TypeName("devolatilisationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DaemDevolatilisationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        DaemDevolatilisationModel(CloudType& owner);

        //- Construct from dictionary
        DaemDevolatilisationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        DaemDevolatilisationModel(const DaemDevolatilisationModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DaemDevolatilisationModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~DaemDevolatilisationModel();


    //- Selector
    static autoPtr<DaemDevolatilisationModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


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
        ) const = 0;

        //- Add to devolatilisation mass
        void addToDaemDevolatilisationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDaemDevolatilisationModel(CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::daemReactingMultiphaseCloudType                   \
        daemReactingMultiphaseCloudType;                                       \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::DaemDevolatilisationModel<daemReactingMultiphaseCloudType>,      \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            DaemDevolatilisationModel<daemReactingMultiphaseCloudType>,        \
            dictionary                                                         \
        );                                                                     \
    }


#define makeDaemDevolatilisationModelType(SS, CloudType)                       \
                                                                               \
    typedef Foam::CloudType::daemReactingMultiphaseCloudType                   \
        daemReactingMultiphaseCloudType;                                       \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<daemReactingMultiphaseCloudType>, 0);                        \
                                                                               \
    Foam::DaemDevolatilisationModel<daemReactingMultiphaseCloudType>::         \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<daemReactingMultiphaseCloudType>>                            \
        add##SS##CloudType##daemReactingMultiphaseCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DaemDevolatilisationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
