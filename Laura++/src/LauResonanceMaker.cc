
/*
Copyright 2004 University of Warwick

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*
Laura++ package authors:
John Back
Paul Harrison
Thomas Latham
*/

/*! \file LauResonanceMaker.cc
    \brief File containing implementation of LauResonanceMaker class.
*/

#include "LauResonanceMaker.hh"

#include "LauAbsResonance.hh"
#include "LauBelleNR.hh"
#include "LauBelleSymNR.hh"
#include "LauBreitWignerRes.hh"
#include "LauDabbaRes.hh"
#include "LauDaughters.hh"
#include "LauEFKLLMRes.hh"
#include "LauFlatNR.hh"
#include "LauFlatteRes.hh"
#include "LauGaussIncohRes.hh"
#include "LauGounarisSakuraiRes.hh"
#include "LauKappaRes.hh"
#include "LauLASSBWRes.hh"
#include "LauLASSNRRes.hh"
#include "LauLASSRes.hh"
#include "LauModIndPartWaveMagPhase.hh"
#include "LauModIndPartWaveRealImag.hh"
#include "LauNRAmplitude.hh"
#include "LauPolNR.hh"
#include "LauPolarFormFactorNR.hh"
#include "LauPolarFormFactorSymNR.hh"
#include "LauPoleRes.hh"
#include "LauRelBreitWignerRes.hh"
#include "LauRescattering2Res.hh"
#include "LauRescatteringRes.hh"
#include "LauResonanceInfo.hh"
#include "LauRhoOmegaMix.hh"
#include "LauSigmaRes.hh"

#include <iostream>

LauResonanceMaker* LauResonanceMaker::resonanceMaker_ = 0;

LauResonanceMaker::LauResonanceMaker() :
    nResDefMax_( 0 ),
    bwBarrierType_( LauBlattWeisskopfFactor::BWPrimeBarrier ),
    bwRestFrame_( LauBlattWeisskopfFactor::ResonanceFrame ),
    spinFormalism_( LauAbsResonance::Zemach_P ),
    summaryPrinted_( kFALSE )
{
    this->createResonanceVector();
    this->setDefaultBWRadius( LauBlattWeisskopfFactor::Parent, 4.0 );
}

LauResonanceMaker::~LauResonanceMaker()
{
    for ( std::vector<LauBlattWeisskopfFactor*>::iterator iter = bwIndepFactors_.begin();
          iter != bwIndepFactors_.end();
          ++iter ) {
        delete *iter;
    }
    bwIndepFactors_.clear();
    for ( BWFactorCategoryMap::iterator iter = bwFactors_.begin(); iter != bwFactors_.end(); ++iter ) {
        delete iter->second.bwFactor_;
    }
    bwFactors_.clear();
}

LauResonanceMaker& LauResonanceMaker::get()
{
    if ( resonanceMaker_ == 0 ) {
        resonanceMaker_ = new LauResonanceMaker();
    }

    return *resonanceMaker_;
}

void LauResonanceMaker::createResonanceVector()
{
    // Function to create all possible resonances that this class supports.
    // Also add in the sigma and kappa - but a special paramterisation is used
    // instead of the PDG "pole mass and width" values.

    std::cout << "INFO in LauResonanceMaker::createResonanceVector : Setting up possible resonance states..."
              << std::endl;

    LauResonanceInfo* neutral( 0 );
    LauResonanceInfo* positve( 0 );
    LauResonanceInfo* negatve( 0 );

    // Define the resonance names and store them in the array
    resInfo_.clear();
    resInfo_.reserve( 100 );
    // rho resonances              name,            mass,     width,    spin,  charge,  default BW category,             BW radius parameter (defaults to 4.0)
    // rho(770)
    neutral =
        new LauResonanceInfo( "rho0(770)", 0.77526, 0.1478, 1, 0, LauBlattWeisskopfFactor::Light, 5.3 );
    positve =
        new LauResonanceInfo( "rho+(770)", 0.77511, 0.1491, 1, 1, LauBlattWeisskopfFactor::Light, 5.3 );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // The following two lines of code are placed here in order to allow the following, rather niche, scenario:
    // The LauRhoOmegaMix code permits (through the use of the optional independentPar argument of LauResonanceInfo::addExtraParameter) the magnitude and phase of the rho/omega mixing to potentially differ between the decay of the parent particle to rho0 X and the parent antiparticle to rho0 Xbar.
    // This can be acheived by using the rho0(770) record in one case and the rho0(770)_COPY record in the other.
    neutral = neutral->createSharedParameterRecord( "rho0(770)_COPY" );
    resInfo_.push_back( neutral );
    // rho(1450)
    neutral = new LauResonanceInfo( "rho0(1450)", 1.465, 0.400, 1, 0, LauBlattWeisskopfFactor::Light );
    positve = new LauResonanceInfo( "rho+(1450)", 1.465, 0.400, 1, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // rho_3(1690)
    neutral =
        new LauResonanceInfo( "rho0_3(1690)", 1.686, 0.186, 3, 0, LauBlattWeisskopfFactor::Light );
    positve =
        new LauResonanceInfo( "rho+_3(1690)", 1.686, 0.186, 3, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // rho(1700)
    neutral = new LauResonanceInfo( "rho0(1700)", 1.720, 0.250, 1, 0, LauBlattWeisskopfFactor::Light );
    positve = new LauResonanceInfo( "rho+(1700)", 1.720, 0.250, 1, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // rho(1900)
    neutral = new LauResonanceInfo( "rho0(1900)", 1.909, 0.130, 1, 0, LauBlattWeisskopfFactor::Light );
    positve = new LauResonanceInfo( "rho+(1900)", 1.909, 0.130, 1, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // rho_3(1990)
    neutral =
        new LauResonanceInfo( "rho0_3(1990)", 1.982, 0.188, 3, 0, LauBlattWeisskopfFactor::Light );
    positve =
        new LauResonanceInfo( "rho+_3(1990)", 1.982, 0.188, 3, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );

    // K* resonances               name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // K*(892)
    neutral =
        new LauResonanceInfo( "K*0(892)", 0.89555, 0.0473, 1, 0, LauBlattWeisskopfFactor::Kstar, 3.0 );
    positve =
        new LauResonanceInfo( "K*+(892)", 0.89166, 0.0508, 1, 1, LauBlattWeisskopfFactor::Kstar, 3.0 );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // K*(1410)
    neutral = new LauResonanceInfo( "K*0(1410)", 1.414, 0.232, 1, 0, LauBlattWeisskopfFactor::Kstar );
    positve = new LauResonanceInfo( "K*+(1410)", 1.414, 0.232, 1, 1, LauBlattWeisskopfFactor::Kstar );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // K*_0(1430)
    neutral = new LauResonanceInfo( "K*0_0(1430)", 1.425, 0.270, 0, 0, LauBlattWeisskopfFactor::Kstar );
    positve = new LauResonanceInfo( "K*+_0(1430)", 1.425, 0.270, 0, 1, LauBlattWeisskopfFactor::Kstar );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // LASS nonresonant model
    neutral = neutral->createSharedParameterRecord( "LASSNR0" );
    positve = positve->createSharedParameterRecord( "LASSNR+" );
    negatve = negatve->createSharedParameterRecord( "LASSNR-" );
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // K*_2(1430)
    neutral =
        new LauResonanceInfo( "K*0_2(1430)", 1.4324, 0.109, 2, 0, LauBlattWeisskopfFactor::Kstar );
    positve =
        new LauResonanceInfo( "K*+_2(1430)", 1.4273, 0.100, 2, 1, LauBlattWeisskopfFactor::Kstar );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // K*(1680)
    neutral = new LauResonanceInfo( "K*0(1680)", 1.718, 0.322, 1, 0, LauBlattWeisskopfFactor::Kstar );
    positve = new LauResonanceInfo( "K*+(1680)", 1.718, 0.322, 1, 1, LauBlattWeisskopfFactor::Kstar );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // K*(1950)
    neutral = new LauResonanceInfo( "K*0_0(1950)", 1.945, 0.201, 0, 0, LauBlattWeisskopfFactor::Kstar );
    positve = new LauResonanceInfo( "K*+_0(1950)", 1.945, 0.201, 0, 1, LauBlattWeisskopfFactor::Kstar );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );

    // phi resonances              name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // phi(1020)
    neutral =
        new LauResonanceInfo( "phi(1020)", 1.019461, 0.004249, 1, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // phi(1680)
    neutral = new LauResonanceInfo( "phi(1680)", 1.680, 0.150, 1, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );

    // f resonances                name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // f_0(980)
    neutral = new LauResonanceInfo( "f_0(980)", 0.990, 0.070, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(1270)
    neutral = new LauResonanceInfo( "f_2(1270)", 1.2755, 0.1867, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_0(1370)
    neutral = new LauResonanceInfo( "f_0(1370)", 1.370, 0.350, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f'_0(1300), from Belle's Kspipi paper
    neutral = new LauResonanceInfo( "f'_0(1300)", 1.449, 0.126, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(1430)
    neutral = new LauResonanceInfo( "f_2(1430)",
                                    1.430,
                                    0.150,
                                    2,
                                    0,
                                    LauBlattWeisskopfFactor::Light );    // PDG width in the range 13 - 150
    resInfo_.push_back( neutral );
    // f_0(1500)
    neutral = new LauResonanceInfo( "f_0(1500)", 1.506, 0.112, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f'_2(1525)
    neutral = new LauResonanceInfo( "f'_2(1525)", 1.5174, 0.086, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(1565)
    neutral = new LauResonanceInfo( "f_2(1565)", 1.542, 0.122, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(1640)
    neutral = new LauResonanceInfo( "f_2(1640)", 1.639, 0.099, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_0(1710)
    neutral = new LauResonanceInfo( "f_0(1710)", 1.704, 0.123, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(1810)
    neutral = new LauResonanceInfo( "f_2(1810)", 1.815, 0.197, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(1910)
    neutral = new LauResonanceInfo( "f_2(1910)", 1.900, 0.167, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(1950)
    neutral = new LauResonanceInfo( "f_2(1950)", 1.936, 0.464, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_2(2010)
    neutral = new LauResonanceInfo( "f_2(2010)", 2.011, 0.202, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_0(2020)
    neutral = new LauResonanceInfo( "f_0(2020)", 1.992, 0.442, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_4(2050)
    neutral = new LauResonanceInfo( "f_4(2050)", 2.018, 0.237, 4, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // f_0(2100)
    neutral = new LauResonanceInfo( "f_0(2100)", 2.086, 0.284, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );

    // omega resonances            name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // omega(782)
    neutral =
        new LauResonanceInfo( "omega(782)", 0.78265, 0.00849, 1, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );

    // a resonances                name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // a_0(980)
    neutral = new LauResonanceInfo( "a0_0(980)", 0.980, 0.092, 0, 0, LauBlattWeisskopfFactor::Light );
    positve = new LauResonanceInfo( "a+_0(980)", 0.980, 0.092, 0, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // a_0(1450)
    neutral = new LauResonanceInfo( "a0_0(1450)", 1.474, 0.265, 0, 0, LauBlattWeisskopfFactor::Light );
    positve = new LauResonanceInfo( "a+_0(1450)", 1.474, 0.265, 0, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // a_2(1320)
    neutral =
        new LauResonanceInfo( "a0_2(1320)", 1.3169, 0.1050, 2, 0, LauBlattWeisskopfFactor::Light );
    positve =
        new LauResonanceInfo( "a+_2(1320)", 1.3169, 0.1050, 2, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );

    // charmonium resonances       name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // chi_c0
    neutral =
        new LauResonanceInfo( "chi_c0", 3.41471, 0.0108, 0, 0, LauBlattWeisskopfFactor::Charmonium );
    resInfo_.push_back( neutral );
    // chi_c1
    neutral =
        new LauResonanceInfo( "chi_c1", 3.51067, 0.00084, 0, 0, LauBlattWeisskopfFactor::Charmonium );
    resInfo_.push_back( neutral );
    // chi_c2
    neutral =
        new LauResonanceInfo( "chi_c2", 3.55617, 0.00197, 2, 0, LauBlattWeisskopfFactor::Charmonium );
    resInfo_.push_back( neutral );
    // psi(3770)
    neutral =
        new LauResonanceInfo( "psi(3770)", 3.7737, 0.0272, 1, 0, LauBlattWeisskopfFactor::Charmonium );
    resInfo_.push_back( neutral );
    // X(3872)
    neutral =
        new LauResonanceInfo( "X(3872)", 3.87169, 0.0012, 1, 0, LauBlattWeisskopfFactor::Charmonium );
    resInfo_.push_back( neutral );
    // chi_c2(2P)
    neutral =
        new LauResonanceInfo( "chi_c2(2P)", 3.9222, 0.0353, 2, 0, LauBlattWeisskopfFactor::Charmonium );
    resInfo_.push_back( neutral );

    // unknown scalars             name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // sigma
    neutral = new LauResonanceInfo( "sigma0", 0.475, 0.550, 0, 0, LauBlattWeisskopfFactor::Light );
    positve = new LauResonanceInfo( "sigma+", 0.475, 0.550, 0, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // kappa
    neutral = new LauResonanceInfo( "kappa0", 0.824, 0.478, 0, 0, LauBlattWeisskopfFactor::Kstar );
    positve = new LauResonanceInfo( "kappa+", 0.824, 0.478, 0, 1, LauBlattWeisskopfFactor::Kstar );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // dabba
    neutral = new LauResonanceInfo( "dabba0", 2.098, 0.520, 0, 0, LauBlattWeisskopfFactor::Charm );
    positve = new LauResonanceInfo( "dabba+", 2.098, 0.520, 0, 1, LauBlattWeisskopfFactor::Charm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );

    // excited charm states        name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // D*
    neutral = new LauResonanceInfo( "D*0", 2.00685, 0.0021, 1, 0, LauBlattWeisskopfFactor::Charm );
    positve = new LauResonanceInfo( "D*+", 2.01026, 83.4e-6, 1, 1, LauBlattWeisskopfFactor::Charm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // D*_0
    neutral = new LauResonanceInfo( "D*0_0", 2.300, 0.274, 0, 0, LauBlattWeisskopfFactor::Charm );
    positve = new LauResonanceInfo( "D*+_0", 2.349, 0.221, 0, 1, LauBlattWeisskopfFactor::Charm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // D*_2
    //AVERAGE--neutral = new LauResonanceInfo("D*0_2",         2.4618,   0.049,    2,     0           );
    neutral = new LauResonanceInfo( "D*0_2", 2.4607, 0.0475, 2, 0, LauBlattWeisskopfFactor::Charm );
    positve = new LauResonanceInfo( "D*+_2", 2.4654, 0.0467, 2, 1, LauBlattWeisskopfFactor::Charm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // D1(2420)
    neutral =
        new LauResonanceInfo( "D0_1(2420)", 2.4208, 0.0317, 1, 0, LauBlattWeisskopfFactor::Charm );
    positve = new LauResonanceInfo( "D+_1(2420)", 2.4232, 0.025, 1, 1, LauBlattWeisskopfFactor::Charm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // D(2600)
    //OLD--neutral = new LauResonanceInfo("D0(2600)",      2.6087,   0.093,    0,     0           );
    //OLD--positve = new LauResonanceInfo("D+(2600)",      2.6213,   0.093,    0,     1           );
    neutral = new LauResonanceInfo( "D0(2600)", 2.623, 0.139, 0, 0, LauBlattWeisskopfFactor::Charm );
    positve = new LauResonanceInfo( "D+(2600)", 2.623, 0.139, 0, 1, LauBlattWeisskopfFactor::Charm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // D(2680)
    neutral =
        new LauResonanceInfo( "D*0_1(2680)", 2.6811, 0.1867, 1, 0, LauBlattWeisskopfFactor::Charm );
    resInfo_.push_back( neutral );
    // D(2760)
    //OLD--    neutral = new LauResonanceInfo("D0(2760)",      2.7633,   0.061,    1,     0           );
    //OLD--    positve = new LauResonanceInfo("D+(2760)",      2.7697,   0.061,    1,     1           );
    neutral = new LauResonanceInfo( "D0(2760)", 2.761, 0.063, 1, 0, LauBlattWeisskopfFactor::Charm );
    positve = new LauResonanceInfo( "D+(2760)", 2.761, 0.063, 1, 1, LauBlattWeisskopfFactor::Charm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    neutral =
        new LauResonanceInfo( "D*0_3(2760)", 2.7755, 0.0953, 3, 0, LauBlattWeisskopfFactor::Charm );
    resInfo_.push_back( neutral );
    // D(2900)
    neutral = new LauResonanceInfo( "D0(3000)", 3.214, 0.186, 0, 0, LauBlattWeisskopfFactor::Charm );
    resInfo_.push_back( neutral );
    // D(3400)
    neutral = new LauResonanceInfo( "D0(3400)", 3.4, 0.15, 0, 0, LauBlattWeisskopfFactor::Charm );
    resInfo_.push_back( neutral );

    // excited strange charm       name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // Ds*
    positve =
        new LauResonanceInfo( "Ds*+", 2.1121, 0.0019, 1, 1, LauBlattWeisskopfFactor::StrangeCharm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // Ds0*(2317)
    positve = new LauResonanceInfo( "Ds*+_0(2317)",
                                    2.3178,
                                    0.0038,
                                    0,
                                    1,
                                    LauBlattWeisskopfFactor::StrangeCharm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // Ds2*(2573)
    positve = new LauResonanceInfo( "Ds*+_2(2573)",
                                    2.5691,
                                    0.0169,
                                    2,
                                    1,
                                    LauBlattWeisskopfFactor::StrangeCharm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // Ds1*(2700)
    positve = new LauResonanceInfo( "Ds*+_1(2700)",
                                    2.7083,
                                    0.120,
                                    1,
                                    1,
                                    LauBlattWeisskopfFactor::StrangeCharm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // Ds1*(2860)
    positve = new LauResonanceInfo( "Ds*+_1(2860)",
                                    2.859,
                                    0.159,
                                    1,
                                    1,
                                    LauBlattWeisskopfFactor::StrangeCharm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // Ds3*(2860)
    positve = new LauResonanceInfo( "Ds*+_3(2860)",
                                    2.860,
                                    0.053,
                                    3,
                                    1,
                                    LauBlattWeisskopfFactor::StrangeCharm );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );

    // excited bottom states       name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // B*
    neutral = new LauResonanceInfo( "B*0", 5.3247, 0.00, 1, 0, LauBlattWeisskopfFactor::Beauty, 6.0 );
    positve = new LauResonanceInfo( "B*+", 5.3247, 0.00, 1, 1, LauBlattWeisskopfFactor::Beauty, 6.0 );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( neutral );
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );

    // excited strange bottom      name,            mass,     width,    spin,  charge,  BW category,    BW radius parameter (defaults to 4.0)
    // Bs*
    neutral =
        new LauResonanceInfo( "Bs*0", 5.4154, 0.00, 1, 0, LauBlattWeisskopfFactor::StrangeBeauty, 6.0 );
    resInfo_.push_back( neutral );

    // nonresonant models          name,            mass,     width,   spin,   charge,  BW category,    BW radius parameter (defaults to 4.0)
    // Phase-space nonresonant model
    neutral = new LauResonanceInfo( "NonReson", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // Theory-based nonresonant model
    neutral = new LauResonanceInfo( "NRModel", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // Belle nonresonant models
    neutral = new LauResonanceInfo( "BelleSymNR", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "BelleNR", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    positve = new LauResonanceInfo( "BelleNR+", 0.0, 0.0, 0, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    neutral = new LauResonanceInfo( "BelleNR_Swave", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    positve = new LauResonanceInfo( "BelleNR_Swave+", 0.0, 0.0, 0, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    neutral = new LauResonanceInfo( "BelleNR_Pwave", 0.0, 0.0, 1, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    positve = new LauResonanceInfo( "BelleNR_Pwave+", 0.0, 0.0, 1, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    neutral = new LauResonanceInfo( "BelleNR_Dwave", 0.0, 0.0, 2, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    positve = new LauResonanceInfo( "BelleNR_Dwave+", 0.0, 0.0, 2, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    neutral = new LauResonanceInfo( "BelleNR_Fwave", 0.0, 0.0, 3, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    positve = new LauResonanceInfo( "BelleNR_Fwave+", 0.0, 0.0, 3, 1, LauBlattWeisskopfFactor::Light );
    negatve = positve->createChargeConjugate();
    resInfo_.push_back( positve );
    resInfo_.push_back( negatve );
    // Taylor expansion nonresonant model
    neutral = new LauResonanceInfo( "NRTaylor", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    // Polynomial nonresonant models
    neutral = new LauResonanceInfo( "PolNR_S0", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "PolNR_S1", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "PolNR_S2", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "PolNR_P0", 0.0, 0.0, 1, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "PolNR_P1", 0.0, 0.0, 1, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "PolNR_P2", 0.0, 0.0, 1, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );

    // Fake resonances for S-Wave splines
    neutral = new LauResonanceInfo( "Spline_S0", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "Spline_S0_Bar", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );

    // Polar Form Factor  nonresonant model
    neutral = new LauResonanceInfo( "PolarFFSymNR", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );
    neutral = new LauResonanceInfo( "PolarFFNR", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );

    // PiPi-KK Inelastic Scattering
    neutral = new LauResonanceInfo( "Rescattering", 0.0, 0.0, 0, 0, LauBlattWeisskopfFactor::Light );
    resInfo_.push_back( neutral );

    nResDefMax_ = resInfo_.size();
}

void LauResonanceMaker::setBWType( const LauBlattWeisskopfFactor::BarrierType bwType )
{
    // Check whether any BW factors have been created and bail out if so
    if ( ! bwIndepFactors_.empty() ) {
        std::cerr << "ERROR in LauResonanceMaker::setBWType : some barrier factors have already been created - cannot change the barrier type now!"
                  << std::endl;
        return;
    }
    for ( BWFactorCategoryMap::const_iterator iter = bwFactors_.begin(); iter != bwFactors_.end();
          ++iter ) {
        if ( iter->second.bwFactor_ != 0 ) {
            std::cerr << "ERROR in LauResonanceMaker::setBWType : some barrier factors have already been created - cannot change the barrier type now!"
                      << std::endl;
            return;
        }
    }

    bwBarrierType_ = bwType;
}

void LauResonanceMaker::setBWBachelorRestFrame( const LauBlattWeisskopfFactor::RestFrame restFrame )
{
    // Check whether any BW factors have been created and bail out if so
    if ( ! bwIndepFactors_.empty() ) {
        std::cerr << "ERROR in LauResonanceMaker::setBWBachelorRestFrame : some barrier factors have already been created - cannot change the rest frame now!"
                  << std::endl;
        return;
    }
    for ( BWFactorCategoryMap::const_iterator iter = bwFactors_.begin(); iter != bwFactors_.end();
          ++iter ) {
        if ( iter->second.bwFactor_ != 0 ) {
            std::cerr << "ERROR in LauResonanceMaker::setBWBachelorRestFrame : some barrier factors have already been created - cannot change the rest frame now!"
                      << std::endl;
            return;
        }
    }

    bwRestFrame_ = restFrame;
}

void LauResonanceMaker::setSpinFormalism( const LauAbsResonance::LauSpinType spinType )
{
    if ( summaryPrinted_ ) {
        std::cerr << "ERROR in LauResonanceMaker::setSpinFormalism : cannot redefine the spin formalism after creating one or more resonances"
                  << std::endl;
        return;
    }
    spinFormalism_ = spinType;
}

void LauResonanceMaker::setDefaultBWRadius(
    const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
    const Double_t bwRadius )
{
    if ( bwCategory == LauBlattWeisskopfFactor::Default ||
         bwCategory == LauBlattWeisskopfFactor::Indep ) {
        std::cerr << "WARNING in LauResonanceMaker::setDefaultBWRadius : cannot set radius values for Default or Indep categories"
                  << std::endl;
        return;
    }

    // Check if we have an information object for this category
    BWFactorCategoryMap::iterator factor_iter = bwFactors_.find( bwCategory );
    if ( factor_iter != bwFactors_.end() ) {
        // If so, we can set the value in the information object
        BlattWeisskopfCategoryInfo& categoryInfo = factor_iter->second;
        categoryInfo.defaultRadius_              = bwRadius;

        // Then we can check if a LauBlattWeisskopfFactor object has been created for this category
        LauBlattWeisskopfFactor* bwFactor = categoryInfo.bwFactor_;
        if ( bwFactor != 0 ) {
            // If it has then we can also set its radius value directly
            LauParameter* radius = bwFactor->getRadiusParameter();
            radius->value( bwRadius );
            radius->initValue( bwRadius );
            radius->genValue( bwRadius );
        }
    } else {
        // If not then we just store the value to be used later
        BlattWeisskopfCategoryInfo& categoryInfo = bwFactors_[bwCategory];
        categoryInfo.bwFactor_                   = 0;
        categoryInfo.defaultRadius_              = bwRadius;
        categoryInfo.radiusFixed_                = kTRUE;
    }
}

void LauResonanceMaker::fixBWRadius( const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
                                     const Bool_t fixRadius )
{
    if ( bwCategory == LauBlattWeisskopfFactor::Default ||
         bwCategory == LauBlattWeisskopfFactor::Indep ) {
        std::cerr << "WARNING in LauResonanceMaker::fixBWRadius : cannot fix/float radius values for Default or Indep categories"
                  << std::endl;
        return;
    }

    // Check if we have an information object for this category
    BWFactorCategoryMap::iterator factor_iter = bwFactors_.find( bwCategory );
    if ( factor_iter != bwFactors_.end() ) {
        // If so, we can set the value in the information object
        BlattWeisskopfCategoryInfo& categoryInfo = factor_iter->second;
        categoryInfo.radiusFixed_                = fixRadius;

        // Then we can check if a LauBlattWeisskopfFactor object has been created for this category
        LauBlattWeisskopfFactor* bwFactor = categoryInfo.bwFactor_;
        if ( bwFactor != 0 ) {
            // If it has then we can also fix/float its radius value directly
            LauParameter* radius = bwFactor->getRadiusParameter();
            radius->fixed( fixRadius );
        }
    } else {
        // If not then we just store the value to be used later
        BlattWeisskopfCategoryInfo& categoryInfo = bwFactors_[bwCategory];
        categoryInfo.bwFactor_                   = 0;
        categoryInfo.defaultRadius_              = -1.0;
        categoryInfo.radiusFixed_                = fixRadius;
    }
}

LauBlattWeisskopfFactor* LauResonanceMaker::getParentBWFactor(
    Int_t resSpin,
    LauBlattWeisskopfFactor::BarrierType barrierType )
{
    LauBlattWeisskopfFactor* bwFactor( 0 );

    // Look up the category in the category information map
    BWFactorCategoryMap::iterator factor_iter = bwFactors_.find( LauBlattWeisskopfFactor::Parent );

    if ( factor_iter == bwFactors_.end() ) {
        // If the category is currently undefined we need to create it
        bwFactor = new LauBlattWeisskopfFactor( resSpin,
                                                4.0,
                                                bwBarrierType_,
                                                bwRestFrame_,
                                                LauBlattWeisskopfFactor::Parent );
        std::cerr << "WARNING in LauResonanceMaker::getParentBWFactor : Default radius 4.0 set for Blatt-Weisskopf factor category: Parent"
                  << std::endl;

        BlattWeisskopfCategoryInfo& categoryInfo = bwFactors_[LauBlattWeisskopfFactor::Parent];
        categoryInfo.bwFactor_                   = bwFactor;
        categoryInfo.defaultRadius_              = bwFactor->getRadiusParameter()->value();
        categoryInfo.radiusFixed_                = kTRUE;
    } else {
        // If it exists, we can check if the factor object has been created
        BlattWeisskopfCategoryInfo& categoryInfo = factor_iter->second;

        if ( categoryInfo.bwFactor_ != 0 ) {
            // If so, simply clone it
            bwFactor = categoryInfo.bwFactor_->createClone( resSpin, barrierType );
        } else {
            // Otherwise we need to create it, using the default value if it has been set
            if ( categoryInfo.defaultRadius_ >= 0.0 ) {
                bwFactor = new LauBlattWeisskopfFactor( resSpin,
                                                        categoryInfo.defaultRadius_,
                                                        bwBarrierType_,
                                                        bwRestFrame_,
                                                        LauBlattWeisskopfFactor::Parent );
            } else {
                bwFactor = new LauBlattWeisskopfFactor( resSpin,
                                                        4.0,
                                                        bwBarrierType_,
                                                        bwRestFrame_,
                                                        LauBlattWeisskopfFactor::Parent );
                std::cerr << "WARNING in LauResonanceMaker::getParentBWFactor : Default radius 4.0 set for Blatt-Weisskopf factor category: Parent"
                          << std::endl;
            }
            categoryInfo.bwFactor_ = bwFactor;

            // Set whether the radius should be fixed/floated
            LauParameter* radius = bwFactor->getRadiusParameter();
            radius->fixed( categoryInfo.radiusFixed_ );
        }
    }

    return bwFactor;
}

LauBlattWeisskopfFactor* LauResonanceMaker::getBWFactor(
    const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
    const LauResonanceInfo* resInfo )
{
    LauBlattWeisskopfFactor* bwFactor( 0 );

    // If this is an independent factor, create it and add it to the list of independent factors, then return it
    if ( bwCategory == LauBlattWeisskopfFactor::Indep ) {
        bwFactor = new LauBlattWeisskopfFactor( *resInfo, bwBarrierType_, bwRestFrame_, bwCategory );
        bwIndepFactors_.push_back( bwFactor );
        return bwFactor;
    }

    // Otherwise, look up the category in the category information map
    BWFactorCategoryMap::iterator factor_iter = bwFactors_.find( bwCategory );

    if ( factor_iter == bwFactors_.end() ) {
        // If the category is currently undefined we need to create it
        bwFactor = new LauBlattWeisskopfFactor( *resInfo, bwBarrierType_, bwRestFrame_, bwCategory );

        BlattWeisskopfCategoryInfo& categoryInfo = bwFactors_[bwCategory];
        categoryInfo.bwFactor_                   = bwFactor;
        categoryInfo.defaultRadius_              = bwFactor->getRadiusParameter()->value();
        categoryInfo.radiusFixed_                = kTRUE;
    } else {
        // If it exists, we can check if the factor object has been created
        BlattWeisskopfCategoryInfo& categoryInfo = factor_iter->second;

        if ( categoryInfo.bwFactor_ != 0 ) {
            // If so, simply clone it
            const UInt_t resSpin = resInfo->getSpin();
            bwFactor             = categoryInfo.bwFactor_->createClone( resSpin,
                                                            categoryInfo.bwFactor_->getBarrierType() );
        } else {
            // Otherwise we need to create it, using the default value if it has been set
            if ( categoryInfo.defaultRadius_ >= 0.0 ) {
                bwFactor = new LauBlattWeisskopfFactor( *resInfo,
                                                        categoryInfo.defaultRadius_,
                                                        bwBarrierType_,
                                                        bwRestFrame_,
                                                        bwCategory );
            } else {
                bwFactor =
                    new LauBlattWeisskopfFactor( *resInfo, bwBarrierType_, bwRestFrame_, bwCategory );
            }
            categoryInfo.bwFactor_ = bwFactor;

            // Set whether the radius should be fixed/floated
            LauParameter* radius = bwFactor->getRadiusParameter();
            radius->fixed( categoryInfo.radiusFixed_ );
        }
    }

    return bwFactor;
}

LauAbsResonance* LauResonanceMaker::getResonance(
    const LauDaughters* daughters,
    const TString& resName,
    const Int_t resPairAmpInt,
    const LauAbsResonance::LauResonanceModel resType,
    const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory )
{
    // Routine to return the appropriate LauAbsResonance object given the resonance
    // name (resName), which daughter is the bachelor track (resPairAmpInt = 1,2 or 3),
    // and the resonance type ("BW" = Breit-Wigner, "Flatte" = Flatte distribution).

    // If this is the first resonance we are making, first print a summary of the formalism
    if ( ! summaryPrinted_ ) {
        std::cout << "INFO in LauResonanceMaker::getResonance : Freezing amplitude formalism:"
                  << std::endl;
        switch ( spinFormalism_ ) {
            case LauAbsResonance::Zemach_P :
                std::cout << "                                        : Spin factors use Zemach spin tensors, with bachelor momentum in resonance rest frame"
                          << std::endl;
                break;
            case LauAbsResonance::Zemach_Pstar :
                std::cout << "                                        : Spin factors use Zemach spin tensors, with bachelor momentum in parent rest frame"
                          << std::endl;
                break;
            case LauAbsResonance::Covariant :
                std::cout << "                                        : Spin factors use Covariant spin tensors"
                          << std::endl;
                break;
            case LauAbsResonance::Legendre :
                std::cout << "                                        : Spin factors are just Legendre polynomials"
                          << std::endl;
                break;
        }
        switch ( bwBarrierType_ ) {
            case LauBlattWeisskopfFactor::BWBarrier :
                std::cout << "                                        : Blatt-Weisskopf barrier factors are the 'non-primed' form"
                          << std::endl;
                break;
            case LauBlattWeisskopfFactor::BWPrimeBarrier :
                std::cout << "                                        : Blatt-Weisskopf barrier factors are the 'primed' form"
                          << std::endl;
                break;
            case LauBlattWeisskopfFactor::ExpBarrier :
                std::cout << "                                        : Blatt-Weisskopf barrier factors are the exponential form"
                          << std::endl;
                break;
        }
        switch ( bwRestFrame_ ) {
            case LauBlattWeisskopfFactor::ParentFrame :
                std::cout << "                                        : Blatt-Weisskopf barrier factors use bachelor momentum in parent rest frame"
                          << std::endl;
                break;
            case LauBlattWeisskopfFactor::ResonanceFrame :
                std::cout << "                                        : Blatt-Weisskopf barrier factors use bachelor momentum in resonance rest frame"
                          << std::endl;
                break;
            case LauBlattWeisskopfFactor::Covariant :
                std::cout << "                                        : Blatt-Weisskopf barrier factors use covariant expression"
                          << std::endl;
                break;
        }

        summaryPrinted_ = kTRUE;
    }

    // Loop over all possible resonance states we have defined in
    // createResonanceVector() until we get a match with the name of the resonance

    LauResonanceInfo* resInfo( 0 );
    for ( std::vector<LauResonanceInfo*>::const_iterator iter = resInfo_.begin();
          iter != resInfo_.end();
          ++iter ) {

        if ( resName == ( *iter )->getName() ) {
            // We have recognised the resonance name.
            std::cout << "INFO in LauResonanceMaker::getResonance : Creating resonance: " << resName
                      << std::endl;

            resInfo = ( *iter );

            // stop looping
            break;
        }
    }

    // if we couldn't find the right resonance then we should return a null pointer
    if ( resInfo == 0 ) {
        std::cout << "ERROR in LauResonanceMaker::getResonance : Unable to locate resonance info for: "
                  << resName << std::endl;
        return 0;
    }

    // Now construct the resonance using the specified type
    LauAbsResonance* theResonance( 0 );
    switch ( resType ) {

        case LauAbsResonance::BW :
            // Simple non-relativistic Breit-Wigner
            std::cout << "                                        : Using simple Breit-Wigner lineshape. "
                      << std::endl;
            theResonance = new LauBreitWignerRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::RelBW : {
            // Relativistic Breit-Wigner with Blatt-Weisskopf factors.
            std::cout << "                                        : Using relativistic Breit-Wigner lineshape. "
                      << std::endl;
            theResonance = new LauRelBreitWignerRes( resInfo, resPairAmpInt, daughters );
            LauBlattWeisskopfFactor::BlattWeisskopfCategory parCategory =
                LauBlattWeisskopfFactor::Parent;
            LauBlattWeisskopfFactor::BlattWeisskopfCategory resCategory = bwCategory;
            if ( bwCategory == LauBlattWeisskopfFactor::Default ) {
                resCategory = resInfo->getBWCategory();
            }
            LauBlattWeisskopfFactor* resBWFactor = this->getBWFactor( resCategory, resInfo );
            LauBlattWeisskopfFactor* parBWFactor = this->getBWFactor( parCategory, resInfo );
            theResonance->setBarrierRadii( resBWFactor, parBWFactor );
            break;
        }

        case LauAbsResonance::GS : {
            // Gounaris-Sakurai function to try and model the rho(770) better
            std::cout << "                                        : Using Gounaris-Sakurai lineshape. "
                      << std::endl;
            theResonance = new LauGounarisSakuraiRes( resInfo, resPairAmpInt, daughters );
            LauBlattWeisskopfFactor::BlattWeisskopfCategory parCategory =
                LauBlattWeisskopfFactor::Parent;
            LauBlattWeisskopfFactor::BlattWeisskopfCategory resCategory = bwCategory;
            if ( bwCategory == LauBlattWeisskopfFactor::Default ) {
                resCategory = resInfo->getBWCategory();
            }
            LauBlattWeisskopfFactor* resBWFactor = this->getBWFactor( resCategory, resInfo );
            LauBlattWeisskopfFactor* parBWFactor = this->getBWFactor( parCategory, resInfo );
            theResonance->setBarrierRadii( resBWFactor, parBWFactor );
            break;
        }

        case LauAbsResonance::Flatte :
            // Flatte distribution - coupled channel Breit-Wigner
            std::cout << "                                        : Using Flatte lineshape. "
                      << std::endl;
            theResonance = new LauFlatteRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::Sigma :
            // Sigma model - should only be used for the pi-pi system
            std::cout << "                                        : Using Sigma lineshape. "
                      << std::endl;
            theResonance = new LauSigmaRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::Kappa :
            // Kappa model - should only be used for the K-pi system
            std::cout << "                                        : Using Kappa lineshape. "
                      << std::endl;
            theResonance = new LauKappaRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::Dabba :
            // Dabba model - should only be used for the D-pi system
            std::cout << "                                        : Using Dabba lineshape. "
                      << std::endl;
            theResonance = new LauDabbaRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::LASS :
            // LASS function to try and model the K-pi S-wave better
            std::cout << "                                        : Using LASS lineshape. "
                      << std::endl;
            theResonance = new LauLASSRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::LASS_BW :
            // LASS function to try and model the K-pi S-wave better
            std::cout << "                                        : Using LASS lineshape (resonant part only). "
                      << std::endl;
            theResonance = new LauLASSBWRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::LASS_NR :
            // LASS function to try and model the K-pi S-wave better
            std::cout << "                                        : Using LASS lineshape (nonresonant part only). "
                      << std::endl;
            theResonance = new LauLASSNRRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::EFKLLM :
            // EFKLLM form-factor description of the K-pi S-wave
            std::cout << "                                        : Using EFKLLM lineshape. "
                      << std::endl;
            theResonance = new LauEFKLLMRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::KMatrix :
            // K-matrix description
            std::cerr << "ERROR in LauResonanceMaker::getResonance : K-matrix type specified, which should be separately handled."
                      << std::endl;
            break;

        case LauAbsResonance::FlatNR :
            // uniform NR amplitude - arguments are there to preserve the interface
            std::cout << "                                        : Using uniform NR lineshape. "
                      << std::endl;
            // we override resPairAmpInt here
            theResonance = new LauFlatNR( resInfo, 0, daughters );
            break;

        case LauAbsResonance::NRModel :
            // NR amplitude model - arguments are there to preserve the interface
            std::cout << "                                        : Using NR-model lineshape. "
                      << std::endl;
            // we override resPairAmpInt here
            theResonance = new LauNRAmplitude( resInfo, 0, daughters );
            break;

        case LauAbsResonance::BelleNR :
        case LauAbsResonance::PowerLawNR :
            // Belle NR amplitude model - arguments are there to preserve the interface
            std::cout << "                                        : Using Belle NR lineshape. "
                      << std::endl;
            theResonance = new LauBelleNR( resInfo, resType, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::BelleSymNR :
        case LauAbsResonance::BelleSymNRNoInter :
        case LauAbsResonance::TaylorNR :
            // Belle NR amplitude model - arguments are there to preserve the interface
            std::cout << "                                        : Using Belle symmetric NR lineshape. "
                      << std::endl;
            theResonance = new LauBelleSymNR( resInfo, resType, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::PolNR :
            // Polynomial NR amplitude model - arguments are there to preserve the interface
            std::cout << "                                        : Using polynomial NR lineshape. "
                      << std::endl;
            theResonance = new LauPolNR( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::Pole :
            // Scalar pole model
            std::cout << "                                        :  Using the scalar Pole lineshape.. "
                      << std::endl;
            theResonance = new LauPoleRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::PolarFFNR :
            // Polar Form Factor NR amplitude model - arguments are there to preserve the interface
            std::cout << "                                        :  Using Polar FormFactor NR lineshape.. "
                      << std::endl;
            theResonance = new LauPolarFormFactorNR( resInfo, resType, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::PolarFFSymNR :
        case LauAbsResonance::PolarFFSymNRNoInter :
            // Polar Form Factor NR amplitude model - arguments are there to preserve the interface
            std::cout << "                                        : Using Polar FormFactor Symetric NR lineshape. "
                      << std::endl;
            theResonance = new LauPolarFormFactorSymNR( resInfo, resType, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::Rescattering :
        case LauAbsResonance::RescatteringNoInter :
            // KKPiPi Inelastic Scattering amplitude - arguments are there to preserve the interface
            std::cout << "                                        : KKPiPi Inelastic Scattering amplitude lineshape. "
                      << std::endl;
            theResonance = new LauRescatteringRes( resInfo, resType, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::Rescattering2 :
            // KKPiPi Inelastic Scattering amplitude - arguments are there to preserve the interface
            std::cout << "                                        : KKPiPi Inelastic Scattering amplitude lineshape. "
                      << std::endl;
            theResonance = new LauRescattering2Res( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::MIPW_MagPhase :
            // Model independent partial wave
            std::cout << "                                        : Using model independent partial wave lineshape (magnitude and phase). "
                      << std::endl;
            theResonance = new LauModIndPartWaveMagPhase( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::MIPW_RealImag :
            // Model independent partial wave
            std::cout << "                                        : Using model independent partial wave lineshape (real and imaginary part). "
                      << std::endl;
            theResonance = new LauModIndPartWaveRealImag( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::GaussIncoh :
            // Incoherent Gaussian
            std::cout << "                                        : Using incoherent Gaussian lineshape. "
                      << std::endl;
            theResonance = new LauGaussIncohRes( resInfo, resPairAmpInt, daughters );
            break;

        case LauAbsResonance::RhoOmegaMix_GS :
        case LauAbsResonance::RhoOmegaMix_RBW :
        case LauAbsResonance::RhoOmegaMix_GS_1 :
        case LauAbsResonance::RhoOmegaMix_RBW_1 :
            // Rho-omega mass mixing model
            std::cout << "                                        : Using rho-omega mass mixing lineshape. "
                      << std::endl;
            theResonance = new LauRhoOmegaMix( resInfo, resType, resPairAmpInt, daughters );
            LauBlattWeisskopfFactor::BlattWeisskopfCategory parCategory =
                LauBlattWeisskopfFactor::Parent;
            LauBlattWeisskopfFactor::BlattWeisskopfCategory resCategory = bwCategory;
            if ( bwCategory == LauBlattWeisskopfFactor::Default ) {
                resCategory = resInfo->getBWCategory();
            }
            LauBlattWeisskopfFactor* resBWFactor = this->getBWFactor( resCategory, resInfo );
            LauBlattWeisskopfFactor* parBWFactor = this->getBWFactor( parCategory, resInfo );
            theResonance->setBarrierRadii( resBWFactor, parBWFactor );
            break;
    }

    // Set the spin formalism choice
    theResonance->setSpinType( spinFormalism_ );

    return theResonance;
}

Int_t LauResonanceMaker::resTypeInt( const TString& name ) const
{
    // Internal function that returns the resonance integer, specified by the
    // order of the resonance vector defined in createResonanceVector(),
    // for a given resonance name.
    Int_t resTypInt( -99 );
    Int_t i( 0 );

    for ( std::vector<LauResonanceInfo*>::const_iterator iter = resInfo_.begin();
          iter != resInfo_.end();
          ++iter ) {

        if ( name.BeginsWith( ( *iter )->getName(), TString::kExact ) == kTRUE ) {
            // We have recognised the resonance from those that are available
            resTypInt = i;
            return resTypInt;
        }
        ++i;
    }

    return resTypInt;
}

void LauResonanceMaker::printAll( std::ostream& stream ) const
{
    for ( std::vector<LauResonanceInfo*>::const_iterator iter = resInfo_.begin();
          iter != resInfo_.end();
          ++iter ) {
        stream << ( **iter ) << std::endl;
    }
}

LauResonanceInfo* LauResonanceMaker::getResInfo( const TString& resName ) const
{
    LauResonanceInfo* resInfo( 0 );
    for ( std::vector<LauResonanceInfo*>::const_iterator iter = resInfo_.begin();
          iter != resInfo_.end();
          ++iter ) {

        if ( resName == ( *iter )->getName() ) {
            // We have recognised the resonance name.
            resInfo = ( *iter );
            // stop looping
            break;
        }
    }

    return resInfo;
}
