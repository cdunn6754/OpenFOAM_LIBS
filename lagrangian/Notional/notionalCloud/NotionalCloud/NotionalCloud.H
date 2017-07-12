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
    Foam::NotionalCloud

Description
    Templated base class for multiphase reacting cloud

    - Adds to reacting cloud
      - multiphase composition
      - devolatilisatsion
      - surface reactions

SourceFiles
    NotionalCloudI.H
    NotionalCloud.C

\*---------------------------------------------------------------------------*/

#ifndef NotionalCloud_H
#define NotionalCloud_H

#include "notionalCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class DevolatilisationModel;

template<class CloudType>
class SurfaceReactionModel;

/*---------------------------------------------------------------------------*\
                   Class NotionalCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NotionalCloud
:
    public CloudType,
    public notionalCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef NotionalCloud<CloudType> notionalCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<NotionalCloud<CloudType>> cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        NotionalCloud(const NotionalCloud&);

        //- Disallow default bitwise assignment
        void operator=(const NotionalCloud&);


protected:

    // Protected data

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Devolatilisation model
            autoPtr
            <
                DevolatilisationModel<NotionalCloud<CloudType>>
            >
            devolatilisationModel_;

            //- Surface reaction model
            autoPtr
            <
                SurfaceReactionModel<NotionalCloud<CloudType>>
            >
            surfaceReactionModel_;


        // Check

            //- Total mass transferred to continuous phase via devolatilisation
            scalar dMassDevolatilisation_;

            //- Total mass transferred to continuous phase via surface
            //  reactions
            scalar dMassSurfaceReaction_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(NotionalCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        NotionalCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );


        //- Copy constructor with new name
        NotionalCloud
        (
            NotionalCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        NotionalCloud
        (
            const fvMesh& mesh,
            const word& name,
            const NotionalCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new NotionalCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new NotionalCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~NotionalCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const NotionalCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();


            // Sub-models

                //- Return const access to devolatilisation model
                inline const DevolatilisationModel
                <
                    NotionalCloud<CloudType>
                >&
                devolatilisation() const;

                //- Return reference to devolatilisation model
                inline DevolatilisationModel
                <
                    NotionalCloud<CloudType>
                >&
                devolatilisation();

                //- Return const access to reacting surface reaction model
                inline const SurfaceReactionModel
                <
                    NotionalCloud<CloudType>
                >&
                surfaceReaction() const;

                //- Return reference to reacting surface reaction model
                inline SurfaceReactionModel
                <
                    NotionalCloud<CloudType>
                >&
                surfaceReaction();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Evolve the cloud
            void evolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();

            //- Write the field data for the cloud
            virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NotionalCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NotionalCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
