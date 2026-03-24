
/*
Copyright 2014 University of Warwick

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

/*! \file LauPolarGammaCPCoeffSet.cc
    \brief File containing implementation of LauPolarGammaCPCoeffSet class.
*/

#include "LauPolarGammaCPCoeffSet.hh"

#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauParameter.hh"
#include "LauPrint.hh"

#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>
#include <vector>

LauParameter* LauPolarGammaCPCoeffSet::gammaGlobal_  = 0;
LauParameter* LauPolarGammaCPCoeffSet::rDGlobal_     = 0;
LauParameter* LauPolarGammaCPCoeffSet::deltaDGlobal_ = 0;

LauPolarGammaCPCoeffSet::LauPolarGammaCPCoeffSet( const TString& compName,
                                                  const DecayType decayType,
                                                  const Double_t x,
                                                  const Double_t y,
                                                  const Double_t rB,
                                                  const Double_t deltaB,
                                                  const Double_t gamma,
                                                  const Double_t rD,
                                                  const Double_t deltaD,
                                                  const Bool_t xFixed,
                                                  const Bool_t yFixed,
                                                  const Bool_t rBFixed,
                                                  const Bool_t deltaBFixed,
                                                  const Bool_t gammaFixed,
                                                  const Bool_t rDFixed,
                                                  const Bool_t deltaDFixed,
                                                  const Bool_t rBSecondStage,
                                                  const Bool_t deltaBSecondStage,
                                                  const Bool_t gammaSecondStage,
                                                  const Bool_t rDSecondStage,
                                                  const Bool_t deltaDSecondStage,
                                                  const Bool_t useGlobalGamma,
                                                  const Bool_t useGlobalADSPars ) :
    LauAbsCoeffSet( compName ),
    decayType_( decayType ),
    x_( 0 ),
    y_( 0 ),
    rB_( 0 ),
    deltaB_( 0 ),
    gamma_( 0 ),
    rD_( 0 ),
    deltaD_( 0 ),
    useGlobalGamma_( useGlobalGamma ),
    useGlobalADSPars_( useGlobalADSPars ),
    nonCPPart_( x, y ),
    cpPart_( 0.0, 0.0 ),
    cpAntiPart_( 0.0, 0.0 ),
    particleCoeff_( 0.0, 0.0 ),
    antiparticleCoeff_( 0.0, 0.0 ),
    acp_( "ACP", 0.0, -1.0, 1.0 )
{
    // All of the possible D decay types need these two parameters
    x_ = new LauParameter( "X", x, minRealImagPart_, maxRealImagPart_, xFixed );
    y_ = new LauParameter( "Y", y, minRealImagPart_, maxRealImagPart_, yFixed );

    // if we're using a global gamma, create it if it doesn't already exist then set gamma_ to point to it
    // otherwise create our individual copy of gamma
    if ( useGlobalGamma_ ) {
        if ( ! gammaGlobal_ ) {
            gammaGlobal_ = new LauParameter( "gamma", gamma, minPhase_, maxPhase_, gammaFixed );
            gamma_       = gammaGlobal_;
        } else {
            gamma_ = gammaGlobal_->createClone();
        }
    } else {
        gamma_ = new LauParameter( "gamma", gamma, minPhase_, maxPhase_, gammaFixed );
    }
    if ( gammaSecondStage && ! gammaFixed ) {
        gamma_->secondStage( kTRUE );
        gamma_->initValue( 0.0 );
    }

    // which of the other parameter we need depends on the D-decay type
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        rB_     = new LauParameter( "rB", rB, minMagnitude_, maxMagnitude_, rBFixed );
        deltaB_ = new LauParameter( "deltaB", deltaB, minPhase_, maxPhase_, deltaBFixed );
        if ( rBSecondStage && ! rBFixed ) {
            rB_->secondStage( kTRUE );
            rB_->initValue( 0.0 );
        }
        if ( deltaBSecondStage && ! deltaBFixed ) {
            deltaB_->secondStage( kTRUE );
            deltaB_->initValue( 0.0 );
        }
    }

    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
         decayType_ == ADS_Favoured_btouOnly ) {
        if ( useGlobalADSPars_ ) {
            if ( ! rDGlobal_ ) {
                rDGlobal_ = new LauParameter( "rD", rD, minMagnitude_, maxMagnitude_, rDFixed );
                deltaDGlobal_ = new LauParameter( "deltaD", deltaD, minPhase_, maxPhase_, deltaDFixed );
                rD_     = rDGlobal_;
                deltaD_ = deltaDGlobal_;
            } else {
                rD_     = rDGlobal_->createClone();
                deltaD_ = deltaDGlobal_->createClone();
            }
        } else {
            rD_     = new LauParameter( "rD", rD, minMagnitude_, maxMagnitude_, rDFixed );
            deltaD_ = new LauParameter( "deltaD", deltaD, minPhase_, maxPhase_, deltaDFixed );
        }
        if ( rDSecondStage && ! rDFixed ) {
            rD_->secondStage( kTRUE );
            rD_->initValue( 0.0 );
        }
        if ( deltaDSecondStage && ! deltaDFixed ) {
            deltaD_->secondStage( kTRUE );
            deltaD_->initValue( 0.0 );
        }
    }
}

LauPolarGammaCPCoeffSet::LauPolarGammaCPCoeffSet( const LauPolarGammaCPCoeffSet& rhs,
                                                  CloneOption cloneOption,
                                                  Double_t constFactor ) :
    LauAbsCoeffSet( rhs.name() ),
    decayType_( rhs.decayType_ ),
    x_( 0 ),
    y_( 0 ),
    rB_( 0 ),
    deltaB_( 0 ),
    gamma_( 0 ),
    rD_( 0 ),
    deltaD_( 0 ),
    useGlobalGamma_( rhs.useGlobalGamma_ ),
    useGlobalADSPars_( rhs.useGlobalADSPars_ ),
    nonCPPart_( rhs.nonCPPart_ ),
    cpPart_( rhs.cpPart_ ),
    cpAntiPart_( rhs.cpAntiPart_ ),
    particleCoeff_( rhs.particleCoeff_ ),
    antiparticleCoeff_( rhs.antiparticleCoeff_ ),
    acp_( rhs.acp_ )
{
    if ( cloneOption == All || cloneOption == TieRealPart ) {
        x_ = rhs.x_->createClone( constFactor );
    } else {
        x_ = new LauParameter( "X",
                               rhs.x_->value(),
                               minRealImagPart_,
                               maxRealImagPart_,
                               rhs.x_->fixed() );
        if ( rhs.x_->blind() ) {
            const LauBlind* blinder = rhs.x_->blinder();
            x_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TieImagPart ) {
        y_ = rhs.y_->createClone( constFactor );
    } else {
        y_ = new LauParameter( "Y",
                               rhs.y_->value(),
                               minRealImagPart_,
                               maxRealImagPart_,
                               rhs.y_->fixed() );
        if ( rhs.y_->blind() ) {
            const LauBlind* blinder = rhs.y_->blinder();
            y_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TieCPPars ) {
        gamma_ = rhs.gamma_->createClone( constFactor );
        if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
             decayType_ == GLW_CPOdd || decayType_ == GLW_CPEven ) {
            rB_     = rhs.rB_->createClone( constFactor );
            deltaB_ = rhs.deltaB_->createClone( constFactor );
        }
        if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
             decayType_ == ADS_Favoured_btouOnly ) {
            rD_     = rhs.rD_->createClone( constFactor );
            deltaD_ = rhs.deltaD_->createClone( constFactor );
        }
    } else {
        if ( useGlobalGamma_ ) {
            gamma_ = gammaGlobal_->createClone();
        } else {
            gamma_ = new LauParameter( "gamma",
                                       rhs.gamma_->value(),
                                       minPhase_,
                                       maxPhase_,
                                       rhs.gamma_->fixed() );
            if ( rhs.gamma_->blind() ) {
                const LauBlind* blinder = rhs.gamma_->blinder();
                gamma_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
            }
            if ( rhs.gamma_->secondStage() && ! rhs.gamma_->fixed() ) {
                gamma_->secondStage( kTRUE );
                gamma_->initValue( 0.0 );
            }
        }
        if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
             decayType_ == GLW_CPOdd || decayType_ == GLW_CPEven ) {
            rB_ = new LauParameter( "rB",
                                    rhs.rB_->value(),
                                    minMagnitude_,
                                    maxMagnitude_,
                                    rhs.rB_->fixed() );
            if ( rhs.rB_->blind() ) {
                const LauBlind* blinder = rhs.rB_->blinder();
                rB_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
            }
            deltaB_ = new LauParameter( "deltaB",
                                        rhs.deltaB_->value(),
                                        minPhase_,
                                        maxPhase_,
                                        rhs.deltaB_->fixed() );
            if ( rhs.deltaB_->blind() ) {
                const LauBlind* blinder = rhs.deltaB_->blinder();
                deltaB_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
            }
            if ( rhs.rB_->secondStage() && ! rhs.rB_->fixed() ) {
                rB_->secondStage( kTRUE );
                rB_->initValue( 0.0 );
            }
            if ( rhs.deltaB_->secondStage() && ! rhs.deltaB_->fixed() ) {
                deltaB_->secondStage( kTRUE );
                deltaB_->initValue( 0.0 );
            }
        }
        if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
             decayType_ == ADS_Favoured_btouOnly ) {
            if ( useGlobalADSPars_ ) {
                rD_     = rDGlobal_->createClone();
                deltaD_ = deltaDGlobal_->createClone();
            } else {
                rD_ = new LauParameter( "rD",
                                        rhs.rD_->value(),
                                        minMagnitude_,
                                        maxMagnitude_,
                                        rhs.rD_->fixed() );
                if ( rhs.rD_->blind() ) {
                    const LauBlind* blinder = rhs.rD_->blinder();
                    rD_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
                }
                deltaD_ = new LauParameter( "deltaD",
                                            rhs.deltaD_->value(),
                                            minPhase_,
                                            maxPhase_,
                                            rhs.deltaD_->fixed() );
                if ( rhs.deltaD_->blind() ) {
                    const LauBlind* blinder = rhs.deltaD_->blinder();
                    deltaD_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
                }
                if ( rhs.rD_->secondStage() && ! rhs.rD_->fixed() ) {
                    rD_->secondStage( kTRUE );
                    rD_->initValue( 0.0 );
                }
                if ( rhs.deltaD_->secondStage() && ! rhs.deltaD_->fixed() ) {
                    deltaD_->secondStage( kTRUE );
                    deltaD_->initValue( 0.0 );
                }
            }
        }
    }
}

void LauPolarGammaCPCoeffSet::adjustName( LauParameter* par, const TString& oldBaseName )
{
    if ( ( par == gamma_ && useGlobalGamma_ ) || ( par == rD_ && useGlobalADSPars_ ) ||
         ( par == deltaD_ && useGlobalADSPars_ ) ) {
        // for global parameters we do not want to adjust their names
        return;
    } else {
        LauAbsCoeffSet::adjustName( par, oldBaseName );
    }
}

std::vector<LauParameter*> LauPolarGammaCPCoeffSet::getParameters()
{
    std::vector<LauParameter*> pars;
    pars.push_back( x_ );
    pars.push_back( y_ );
    if ( ! gamma_->fixed() ) {
        pars.push_back( gamma_ );
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        if ( ! rB_->fixed() ) {
            pars.push_back( rB_ );
        }
        if ( ! deltaB_->fixed() ) {
            pars.push_back( deltaB_ );
        }
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
         decayType_ == ADS_Favoured_btouOnly ) {
        if ( ! rD_->fixed() ) {
            pars.push_back( rD_ );
        }
        if ( ! deltaD_->fixed() ) {
            pars.push_back( deltaD_ );
        }
    }
    return pars;
}

void LauPolarGammaCPCoeffSet::printParValues() const
{
    std::cout << "INFO in LauPolarGammaCPCoeffSet::printParValues : Component \"" << this->name()
              << "\" has ";
    std::cout << "x = " << x_->value() << ",\t";
    std::cout << "y = " << y_->value() << ",\t";
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        std::cout << "rB = " << rB_->value() << ",\t";
        std::cout << "deltaB = " << deltaB_->value() << ",\t";
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
         decayType_ == ADS_Favoured_btouOnly ) {
        std::cout << "rD = " << rD_->value() << ",\t";
        std::cout << "deltaD = " << deltaD_->value() << ",\t";
    }
    std::cout << "gamma = " << gamma_->value() << "." << std::endl;
}

void LauPolarGammaCPCoeffSet::printTableHeading( std::ostream& stream ) const
{
    switch ( decayType_ ) {
        case GLW_CPOdd :
            stream << "\\begin{tabular}{|l|c|c|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $r_B$ & $\\delta_B$ & $\\gamma$ \\\\"
                   << std::endl;
            break;
        case GLW_CPEven :
            stream << "\\begin{tabular}{|l|c|c|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $r_B$ & $\\delta_B$ & $\\gamma$ \\\\"
                   << std::endl;
            break;
        case ADS_Favoured :
            stream << "\\begin{tabular}{|l|c|c|c|c|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $r_B$ & $\\delta_B$ & $r_D$ & $\\delta_D$ & $\\gamma$ \\\\"
                   << std::endl;
            break;
        case ADS_Suppressed :
            stream << "\\begin{tabular}{|l|c|c|c|c|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $r_B$ & $\\delta_B$ & $r_D$ & $\\delta_D$ & $\\gamma$ \\\\"
                   << std::endl;
            break;
        case GLW_CPOdd_btouOnly :
            stream << "\\begin{tabular}{|l|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $\\gamma$ \\\\" << std::endl;
            break;
        case GLW_CPEven_btouOnly :
            stream << "\\begin{tabular}{|l|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $\\gamma$ \\\\" << std::endl;
            break;
        case ADS_Favoured_btouOnly :
            stream << "\\begin{tabular}{|l|c|c|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $r_D$ & $\\delta_D$ & $\\gamma$ \\\\"
                   << std::endl;
            break;
        case ADS_Suppressed_btouOnly :
            stream << "\\begin{tabular}{|l|c|c|c|}" << std::endl;
            stream << "\\hline" << std::endl;
            stream << "Component & Real Part & Imaginary Part & $\\gamma$ \\\\" << std::endl;
            break;
    }
    stream << "\\hline" << std::endl;
}

void LauPolarGammaCPCoeffSet::printTableRow( std::ostream& stream ) const
{
    LauPrint print;
    TString resName = this->name();
    resName         = resName.ReplaceAll( "_", "\\_" );
    stream << resName << "  &  $";
    print.printFormat( stream, x_->value() );
    stream << " \\pm ";
    print.printFormat( stream, x_->error() );
    stream << "$  &  $";
    print.printFormat( stream, y_->value() );
    stream << " \\pm ";
    print.printFormat( stream, y_->error() );
    stream << "$  &  $";
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        print.printFormat( stream, rB_->value() );
        stream << " \\pm ";
        print.printFormat( stream, rB_->error() );
        stream << "$  &  $";
        print.printFormat( stream, deltaB_->value() );
        stream << " \\pm ";
        print.printFormat( stream, deltaB_->error() );
        stream << "$  &  $";
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
         decayType_ == ADS_Favoured_btouOnly ) {
        print.printFormat( stream, rD_->value() );
        stream << " \\pm ";
        print.printFormat( stream, rD_->error() );
        stream << "$  &  $";
        print.printFormat( stream, deltaD_->value() );
        stream << " \\pm ";
        print.printFormat( stream, deltaD_->error() );
        stream << "$  &  $";
    }
    print.printFormat( stream, gamma_->value() );
    stream << " \\pm ";
    print.printFormat( stream, gamma_->error() );
    stream << "$ \\\\" << std::endl;
}

void LauPolarGammaCPCoeffSet::randomiseInitValues()
{
    if ( x_->fixed() == kFALSE ) {
        // Choose a value for "X" between -3.0 and 3.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 6.0 - 3.0;
        x_->initValue( value );
        x_->value( value );
    }
    if ( y_->fixed() == kFALSE ) {
        // Choose a value for "Y" between -3.0 and 3.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 6.0 - 3.0;
        y_->initValue( value );
        y_->value( value );
    }
    if ( gamma_->fixed() == kFALSE && gamma_->secondStage() == kFALSE ) {
        // Choose a value for "gamma" between +-pi
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                         LauConstants::pi;
        gamma_->initValue( value );
        gamma_->value( value );
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        if ( rB_->fixed() == kFALSE && rB_->secondStage() == kFALSE ) {
            // Choose a value for "rB" between 0.0 and 2.0
            Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 2.0;
            rB_->initValue( value );
            rB_->value( value );
        }
        if ( deltaB_->fixed() == kFALSE && deltaB_->secondStage() == kFALSE ) {
            // Choose a value for "deltaB" between +- pi
            Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                             LauConstants::pi;
            deltaB_->initValue( value );
            deltaB_->value( value );
        }
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
         decayType_ == ADS_Favoured_btouOnly ) {
        if ( rD_->fixed() == kFALSE && rD_->secondStage() == kFALSE ) {
            // Choose a value for "rD" between 0.0 and 2.0
            Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 2.0;
            rD_->initValue( value );
            rD_->value( value );
        }
        if ( deltaD_->fixed() == kFALSE && deltaD_->secondStage() == kFALSE ) {
            // Choose a value for "deltaD" between +- pi
            Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                             LauConstants::pi;
            deltaD_->initValue( value );
            deltaD_->value( value );
        }
    }
}

void LauPolarGammaCPCoeffSet::finaliseValues()
{
    // retrieve the current values from the parameters
    Double_t gammaVal  = gamma_->value();
    Double_t rBVal     = 0.0;
    Double_t deltaBVal = 0.0;
    Double_t genDeltaB = 0.0;
    Double_t rDVal     = 0.0;
    Double_t deltaDVal = 0.0;
    Double_t genDeltaD = 0.0;
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        rBVal     = rB_->value();
        deltaBVal = deltaB_->value();
        genDeltaB = deltaB_->genValue();
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
         decayType_ == ADS_Favoured_btouOnly ) {
        rDVal     = rD_->value();
        deltaDVal = deltaD_->value();
        genDeltaD = deltaD_->genValue();
    }

    // Check whether we have a negative magnitude.
    // If so make it positive and add pi to the phases.
    if ( rBVal < 0.0 ) {
        rBVal     *= -1.0;
        deltaBVal += LauConstants::pi;
    }
    if ( rDVal < 0.0 ) {
        rDVal     *= -1.0;
        deltaDVal += LauConstants::pi;
    }

    // Check now whether the phases lie in the right range (-pi to pi).
    Bool_t deltaBWithinRange( kFALSE );
    Bool_t deltaDWithinRange( kFALSE );
    Bool_t gammaWithinRange( kFALSE );
    while ( deltaBWithinRange == kFALSE ) {
        if ( deltaBVal > -LauConstants::pi && deltaBVal <= LauConstants::pi ) {
            deltaBWithinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( deltaBVal > LauConstants::pi ) {
                deltaBVal -= LauConstants::twoPi;
            } else if ( deltaBVal <= -LauConstants::pi ) {
                deltaBVal += LauConstants::twoPi;
            }
        }
    }

    while ( deltaDWithinRange == kFALSE ) {
        if ( deltaDVal > -LauConstants::pi && deltaDVal <= LauConstants::pi ) {
            deltaDWithinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( deltaDVal > LauConstants::pi ) {
                deltaDVal -= LauConstants::twoPi;
            } else if ( deltaDVal <= -LauConstants::pi ) {
                deltaDVal += LauConstants::twoPi;
            }
        }
    }

    while ( gammaWithinRange == kFALSE ) {
        if ( gammaVal > -LauConstants::pi && gammaVal <= LauConstants::pi ) {
            gammaWithinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( gammaVal > LauConstants::pi ) {
                gammaVal -= LauConstants::twoPi;
            } else if ( gammaVal <= -LauConstants::pi ) {
                gammaVal += LauConstants::twoPi;
            }
        }
    }

    // To resolve the two-fold ambiguity in gamma and deltaB we require gamma to be in the range 0-pi
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        if ( gammaVal < 0.0 ) {
            if ( deltaBVal <= 0.0 ) {
                gammaVal  += LauConstants::pi;
                deltaBVal += LauConstants::pi;
            } else {
                gammaVal  += LauConstants::pi;
                deltaBVal -= LauConstants::pi;
            }
        }
    }

    // A further problem can occur when the generated phase is close to -pi or pi.
    // The phase can wrap over to the other end of the scale -
    // this leads to artificially large pulls so we wrap it back.
    Double_t diff = deltaBVal - genDeltaB;
    if ( diff > LauConstants::pi ) {
        deltaBVal -= LauConstants::twoPi;
    } else if ( diff < -LauConstants::pi ) {
        deltaBVal += LauConstants::twoPi;
    }

    diff = deltaDVal - genDeltaD;
    if ( diff > LauConstants::pi ) {
        deltaDVal -= LauConstants::twoPi;
    } else if ( diff < -LauConstants::pi ) {
        deltaDVal += LauConstants::twoPi;
    }

    // finally store the new values in the parameters
    // and update the pulls
    gamma_->value( gammaVal );
    gamma_->updatePull();
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed || decayType_ == GLW_CPOdd ||
         decayType_ == GLW_CPEven ) {
        rB_->value( rBVal );
        rB_->updatePull();
        deltaB_->value( deltaBVal );
        deltaB_->updatePull();
    }
    if ( decayType_ == ADS_Favoured || decayType_ == ADS_Suppressed ||
         decayType_ == ADS_Favoured_btouOnly ) {
        rD_->value( rDVal );
        rD_->updatePull();
        deltaD_->value( deltaDVal );
        deltaD_->updatePull();
    }
}

const LauComplex& LauPolarGammaCPCoeffSet::particleCoeff()
{
    this->updateAmplitudes();
    return particleCoeff_;
}

const LauComplex& LauPolarGammaCPCoeffSet::antiparticleCoeff()
{
    this->updateAmplitudes();
    return antiparticleCoeff_;
}

void LauPolarGammaCPCoeffSet::updateAmplitudes()
{
    nonCPPart_.setRealImagPart( x_->unblindValue(), y_->unblindValue() );

    const Double_t gammaVal = gamma_->unblindValue();

    switch ( decayType_ ) {

        case GLW_CPOdd : {
            const Double_t rBVal     = rB_->unblindValue();
            const Double_t deltaBVal = deltaB_->unblindValue();
            cpPart_.setRealImagPart( 1.0 - rBVal * TMath::Cos( deltaBVal + gammaVal ),
                                     -rBVal * TMath::Sin( deltaBVal + gammaVal ) );
            cpAntiPart_.setRealImagPart( 1.0 - rBVal * TMath::Cos( deltaBVal - gammaVal ),
                                         -rBVal * TMath::Sin( deltaBVal - gammaVal ) );
            break;
        }

        case GLW_CPEven : {
            const Double_t rBVal     = rB_->unblindValue();
            const Double_t deltaBVal = deltaB_->unblindValue();
            cpPart_.setRealImagPart( 1.0 + rBVal * TMath::Cos( deltaBVal + gammaVal ),
                                     rBVal * TMath::Sin( deltaBVal + gammaVal ) );
            cpAntiPart_.setRealImagPart( 1.0 + rBVal * TMath::Cos( deltaBVal - gammaVal ),
                                         rBVal * TMath::Sin( deltaBVal - gammaVal ) );
            break;
        }

        case ADS_Favoured : {
            const Double_t rBVal     = rB_->unblindValue();
            const Double_t deltaBVal = deltaB_->unblindValue();
            const Double_t rDVal     = rD_->unblindValue();
            const Double_t deltaDVal = deltaD_->unblindValue();
            cpPart_.setRealImagPart( 1.0 + rBVal * rDVal *
                                               TMath::Cos( deltaBVal - deltaDVal + gammaVal ),
                                     rBVal * rDVal * TMath::Sin( deltaBVal - deltaDVal + gammaVal ) );
            cpAntiPart_.setRealImagPart(
                1.0 + rBVal * rDVal * TMath::Cos( deltaBVal - deltaDVal - gammaVal ),
                rBVal * rDVal * TMath::Sin( deltaBVal - deltaDVal - gammaVal ) );
            break;
        }

        case ADS_Suppressed : {
            const Double_t rBVal     = rB_->unblindValue();
            const Double_t deltaBVal = deltaB_->unblindValue();
            const Double_t rDVal     = rD_->unblindValue();
            const Double_t deltaDVal = deltaD_->unblindValue();
            cpPart_.setRealImagPart(
                rDVal * TMath::Cos( -deltaDVal ) + rBVal * TMath::Cos( deltaBVal + gammaVal ),
                rDVal * TMath::Sin( -deltaDVal ) + rBVal * TMath::Sin( deltaBVal + gammaVal ) );
            cpAntiPart_.setRealImagPart(
                rDVal * TMath::Cos( -deltaDVal ) + rBVal * TMath::Cos( deltaBVal - gammaVal ),
                rDVal * TMath::Sin( -deltaDVal ) + rBVal * TMath::Sin( deltaBVal - gammaVal ) );
            break;
        }

        case GLW_CPOdd_btouOnly :
            nonCPPart_.rescale( -1.0 );
            cpPart_.setRealImagPart( 1.0 * TMath::Cos( gammaVal ), 1.0 * TMath::Sin( gammaVal ) );
            cpAntiPart_.setRealImagPart( 1.0 * TMath::Cos( -gammaVal ),
                                         1.0 * TMath::Sin( -gammaVal ) );
            break;

        case GLW_CPEven_btouOnly :
            cpPart_.setRealImagPart( 1.0 * TMath::Cos( gammaVal ), 1.0 * TMath::Sin( gammaVal ) );
            cpAntiPart_.setRealImagPart( 1.0 * TMath::Cos( -gammaVal ),
                                         1.0 * TMath::Sin( -gammaVal ) );
            break;

        case ADS_Favoured_btouOnly : {
            const Double_t rDVal     = rD_->unblindValue();
            const Double_t deltaDVal = deltaD_->unblindValue();
            cpPart_.setRealImagPart( rDVal * TMath::Cos( -deltaDVal + gammaVal ),
                                     rDVal * TMath::Sin( -deltaDVal + gammaVal ) );
            cpAntiPart_.setRealImagPart( rDVal * TMath::Cos( -deltaDVal - gammaVal ),
                                         rDVal * TMath::Sin( -deltaDVal - gammaVal ) );
            break;
        }

        case ADS_Suppressed_btouOnly :
            cpPart_.setRealImagPart( 1.0 * TMath::Cos( gammaVal ), 1.0 * TMath::Sin( gammaVal ) );
            cpAntiPart_.setRealImagPart( 1.0 * TMath::Cos( -gammaVal ),
                                         1.0 * TMath::Sin( -gammaVal ) );
            break;
    }

    particleCoeff_     = nonCPPart_ * cpPart_;
    antiparticleCoeff_ = nonCPPart_ * cpAntiPart_;
}

void LauPolarGammaCPCoeffSet::setCoeffValues( const LauComplex&, const LauComplex&, Bool_t )
{
    std::cerr << "ERROR in LauPolarGammaCPCoeffSet::setCoeffValues : Method not supported by this class - too many parameters"
              << std::endl;
}

LauParameter LauPolarGammaCPCoeffSet::acp()
{
    // set the name
    TString parName( this->baseName() );
    parName += "_ACP";
    acp_.name( parName );

    // work out the ACP value
    LauComplex nonCPPart( x_->value(), y_->value() );
    LauComplex cpPart;
    LauComplex cpAntiPart;

    const Double_t gammaVal = gamma_->value();

    switch ( decayType_ ) {

        case GLW_CPOdd : {
            const Double_t rBVal     = rB_->value();
            const Double_t deltaBVal = deltaB_->value();
            cpPart.setRealImagPart( 1.0 - rBVal * TMath::Cos( deltaBVal + gammaVal ),
                                    -rBVal * TMath::Sin( deltaBVal + gammaVal ) );
            cpAntiPart.setRealImagPart( 1.0 - rBVal * TMath::Cos( deltaBVal - gammaVal ),
                                        -rBVal * TMath::Sin( deltaBVal - gammaVal ) );
            break;
        }

        case GLW_CPEven : {
            const Double_t rBVal     = rB_->value();
            const Double_t deltaBVal = deltaB_->value();
            cpPart.setRealImagPart( 1.0 + rBVal * TMath::Cos( deltaBVal + gammaVal ),
                                    rBVal * TMath::Sin( deltaBVal + gammaVal ) );
            cpAntiPart.setRealImagPart( 1.0 + rBVal * TMath::Cos( deltaBVal - gammaVal ),
                                        rBVal * TMath::Sin( deltaBVal - gammaVal ) );
            break;
        }

        case ADS_Favoured : {
            const Double_t rBVal     = rB_->value();
            const Double_t deltaBVal = deltaB_->value();
            const Double_t rDVal     = rD_->value();
            const Double_t deltaDVal = deltaD_->value();
            cpPart.setRealImagPart( 1.0 + rBVal * rDVal *
                                              TMath::Cos( deltaBVal - deltaDVal + gammaVal ),
                                    rBVal * rDVal * TMath::Sin( deltaBVal - deltaDVal + gammaVal ) );
            cpAntiPart.setRealImagPart(
                1.0 + rBVal * rDVal * TMath::Cos( deltaBVal - deltaDVal - gammaVal ),
                rBVal * rDVal * TMath::Sin( deltaBVal - deltaDVal - gammaVal ) );
            break;
        }

        case ADS_Suppressed : {
            const Double_t rBVal     = rB_->value();
            const Double_t deltaBVal = deltaB_->value();
            const Double_t rDVal     = rD_->value();
            const Double_t deltaDVal = deltaD_->value();
            cpPart.setRealImagPart(
                rDVal * TMath::Cos( -deltaDVal ) + rBVal * TMath::Cos( deltaBVal + gammaVal ),
                rDVal * TMath::Sin( -deltaDVal ) + rBVal * TMath::Sin( deltaBVal + gammaVal ) );
            cpAntiPart.setRealImagPart(
                rDVal * TMath::Cos( -deltaDVal ) + rBVal * TMath::Cos( deltaBVal - gammaVal ),
                rDVal * TMath::Sin( -deltaDVal ) + rBVal * TMath::Sin( deltaBVal - gammaVal ) );
            break;
        }

        case GLW_CPOdd_btouOnly :
            nonCPPart.rescale( -1.0 );
            cpPart.setRealImagPart( 1.0 * TMath::Cos( gammaVal ), 1.0 * TMath::Sin( gammaVal ) );
            cpAntiPart.setRealImagPart( 1.0 * TMath::Cos( -gammaVal ), 1.0 * TMath::Sin( -gammaVal ) );
            break;

        case GLW_CPEven_btouOnly :
            cpPart.setRealImagPart( 1.0 * TMath::Cos( gammaVal ), 1.0 * TMath::Sin( gammaVal ) );
            cpAntiPart.setRealImagPart( 1.0 * TMath::Cos( -gammaVal ), 1.0 * TMath::Sin( -gammaVal ) );
            break;

        case ADS_Favoured_btouOnly : {
            const Double_t rDVal     = rD_->value();
            const Double_t deltaDVal = deltaD_->value();
            cpPart.setRealImagPart( rDVal * TMath::Cos( -deltaDVal + gammaVal ),
                                    rDVal * TMath::Sin( -deltaDVal + gammaVal ) );
            cpAntiPart.setRealImagPart( rDVal * TMath::Cos( -deltaDVal - gammaVal ),
                                        rDVal * TMath::Sin( -deltaDVal - gammaVal ) );
            break;
        }

        case ADS_Suppressed_btouOnly :
            cpPart.setRealImagPart( 1.0 * TMath::Cos( gammaVal ), 1.0 * TMath::Sin( gammaVal ) );
            cpAntiPart.setRealImagPart( 1.0 * TMath::Cos( -gammaVal ), 1.0 * TMath::Sin( -gammaVal ) );
            break;
    }

    const LauComplex partCoeff = nonCPPart * cpPart;
    const LauComplex antiCoeff = nonCPPart * cpAntiPart;

    const Double_t numer = antiCoeff.abs2() - partCoeff.abs2();
    const Double_t denom = antiCoeff.abs2() + partCoeff.abs2();
    const Double_t value = numer / denom;

    // is it fixed?
    const Bool_t fixed = gamma_->fixed();
    acp_.fixed( fixed );

    // we can't work out the error without the covariance matrix
    const Double_t error( 0.0 );

    // set the value and error
    acp_.valueAndErrors( value, error );

    return acp_;
}

LauAbsCoeffSet* LauPolarGammaCPCoeffSet::createClone( const TString& newName,
                                                      CloneOption cloneOption,
                                                      Double_t constFactor )
{
    LauAbsCoeffSet* clone( 0 );
    if ( cloneOption == All || cloneOption == TieRealPart || cloneOption == TieImagPart ||
         cloneOption == TieCPPars ) {
        clone = new LauPolarGammaCPCoeffSet( *this, cloneOption, constFactor );
        clone->name( newName );
    } else {
        std::cerr << "ERROR in LauPolarGammaCPCoeffSet::createClone : Invalid clone option"
                  << std::endl;
    }
    return clone;
}
