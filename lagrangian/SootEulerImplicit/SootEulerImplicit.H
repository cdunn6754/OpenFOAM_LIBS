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

Description

\*---------------------------------------------------------------------------*/

//using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifndef EulerImplicitSystem_H
#define EulerImplicitSystem_H

#include <math.h>
#include "scalar.H"
#include "scalarField.H"
#include "dimensionSet.H"
#include "psiReactionThermo.H"
#include "fvScalarMatrix.H"
#include "label.H"
#include "tmp.H"

namespace Foam
{

class EulerImplicitSystem
{

    // true constants
    const scalar rho_s = 2000.; // soot density [kg/m^3]
    const scalar Na = 6.232e26; //Avagadros number [particles/kMol]
    const scalar sigma = 1.38e-23; // boltzman [J/K]
    const scalar Ca = 9.; //agglomeration constant [-]
    const scalar n_c = 60.; // minimum particle number for agglomeration [particles]?
    const scalar pi = 3.1415;

    // make the proper dimension sets for the sources
    dimensionSet Y_source_dims; //[kg/(m^3 * s)]
    dimensionSet N_source_dims; // [1/(m^3 * s)]

    // molecular weights
    const scalar MW_soot; // Molecular weight of soot [kg/kmol]
    const scalar MW_c2h2;
    const scalar MW_oh;
    const scalar MW_o2;

    // Thermo fields
    const scalarField& rho;
    const scalarField& T;
    const scalarField& Y_c2h2;
    const scalarField& Y_oh;
    const scalarField& Y_o2;
    
    // Soot fields
    const scalarField& Y_s;
    const scalarField& N_s;

    // Soot source fields for main Yeqn
    scalarField Y_source;
    scalarField N_source;
    

    // Constants are determined from examination of 
    // soot equation source equations while
    // only considering Y_s and N_s to be varible

    // Particle surface area constant
    scalar C_A;

    // Constants for the equations form of 'C_equation_term'
    scalar C_1_1;
    scalar C_1_2;
    scalar C_1_3;
    scalar C_1_4;
    
    scalar C_2_1;
    scalar C_2_2;


    // Set the constants based on new time step species/rho/T
    void setConstants(const label& cellNumber);

    label nEqns() const;

    // The value of the derivatives given Y and N
    // I.E. dY/dt = f(Y,N) with f taken from Kronenburg 
    // CMC paper
    void derivatives
    (
        const scalar& Y, 
        const scalar& N,
        scalarField& derivative
    );

    // Evaluate the newton equations for the euler implicit form of the 
    // two source equations i.e. 'Y^i - Y^(i-1) - dt(dY/dt^i) = 0
    // for a single cell.
    // Assumes
    // Y_i and N_i hold the values at the integration start time.
    // 
    // Y,N - values for Y_soot and N_soot (guessed values most often)
    // dt - Euler implicit integration time step
    // fY, fN - values of the newton method equations evaluated at Y and N
    void newtonEquations
    (
        const scalar& Y,
        const scalar& N,
        const scalar& Y_i,
        const scalar& N_i,
        const scalar dt,
        scalar& fY,
        scalar& fN
    )const ;

    // Find the Jacobian of the newton method system
    // i.e. the derivatives of the functions from
    // newtonEquations
    // again this peforms the calculation for a single cell.
    void jacobian
    (
        const scalar& Y,
        const scalar& N,
        const scalar dt,
        scalarSquareMatrix& J,
        scalarSquareMatrix& invJ
    ) const;

    void invertJacobian
    (
        const scalarSquareMatrix& J,
        scalarSquareMatrix& invJ
    ) const;

    void explicitEuler
    (
        scalar& Y_i,
        scalar& N_i,
        const scalar& dt
    );


    // Integrate over timestep dt.
    // Once Y^(t+dt),N ^(t + dt) are found use them to find the 
    // source terms. Update sources.
    void correct
    (
        const scalar dt,
        const scalar relTol = 1e-6
    );    


public: 

    // Constructor
    EulerImplicitSystem
    (
        const psiReactionThermo& thermo, 
        const basicSpecieMixture& composition,
        const volScalarField& Ns,
        const volScalarField& rho
    );


    // Return the source for the Ysoot equation
    tmp<fvScalarMatrix> sourceY
    (
        const volScalarField& Y_soot
    );

    
    // Return the source for the N_s equation
    tmp<fvScalarMatrix> sourceN
    (
        const volScalarField& N_soot
    );

    // Main function for calculating the sources.
    // Manages substeps of the main time step dt
    void updateSources
    (
        const scalar dt,
        const label subSteps = 3, 
        const scalar relTol = 1e-6
    );


}; // end of class EulerImlicitSystem

} // End namepace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
