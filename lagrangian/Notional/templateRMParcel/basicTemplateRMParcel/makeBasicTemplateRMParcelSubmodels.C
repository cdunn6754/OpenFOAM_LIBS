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

#include "basicTemplateRMCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"
#include "makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
makeParcelCloudFunctionObjects(basicTemplateRMCloud);

// Kinematic sub-models
makeThermoParcelForces(basicTemplateRMCloud);
makeParcelDispersionModels(basicTemplateRMCloud);
makeReactingMultiphaseParcelInjectionModels(basicTemplateRMCloud);
makeParcelPatchInteractionModels(basicTemplateRMCloud);
makeReactingMultiphaseParcelStochasticCollisionModels
(
    basicTemplateRMCloud
);
makeReactingParcelSurfaceFilmModels(basicTemplateRMCloud);

// Thermo sub-models
makeParcelHeatTransferModels(basicTemplateRMCloud);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels
(
    basicTemplateRMCloud
);
makeReactingParcelPhaseChangeModels(basicTemplateRMCloud);

// Reacting multiphase sub-models
makeReactingMultiphaseParcelDevolatilisationModels
(
    basicTemplateRMCloud
);
makeReactingMultiphaseParcelSurfaceReactionModels
(
    basicTemplateRMCloud
);


// Turbulent stuff 07-13-17:
#include "basicTemplateRMCloud.H"

#include "makeParcelTurbulenceDispersionModels.H"
#include "makeThermoParcelTurbulenceForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeThermoParcelTurbulenceForces(basicTemplateRMCloud);
    makeParcelTurbulenceDispersionModels(basicTemplateRMCloud);
}


// Coal parcel (just surface reaction) stuff:
#include "makeCoalParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeCoalParcelSurfaceReactionModels(basicTemplateRMCloud);
}


// ************************************************************************* //
