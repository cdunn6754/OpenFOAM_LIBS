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



        // References to the cloud sub-models

        // Check


    // Protected Member Functions



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


            // Sub-models


        // Cloud evolution functions


            //- Evolve the cloud
            void evolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O


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
