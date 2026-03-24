
/*
Copyright 2008 University of Warwick

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

/*! \file LauKMatrixPropagator.cc
    \brief File containing implementation of LauKMatrixPropagator class.
*/

#include "LauKMatrixPropagator.hh"

#include "LauComplex.hh"
#include "LauKinematics.hh"
#include "LauTextFileParser.hh"

#include "TMath.h"
#include "TSystem.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

LauKMatrixPropagator::LauKMatrixPropagator( const TString& name,
                                            const TString& paramFile,
                                            Int_t resPairAmpInt,
                                            Int_t nChannels,
                                            Int_t nPoles,
                                            Int_t rowIndex ) :
    name_( name ),
    paramFileName_( paramFile ),
    resPairAmpInt_( resPairAmpInt ),
    index_( rowIndex - 1 ),
    nChannels_( nChannels ),
    nPoles_( nPoles )
{
    // Constructor

    // Check that the index is OK
    if ( index_ < 0 || index_ >= nChannels_ ) {
        std::cerr << "ERROR in LauKMatrixPropagator constructor. The rowIndex, which is set to "
                  << rowIndex << ", must be between 1 and the number of channels " << nChannels_
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    this->setParameters( paramFile );
}

LauKMatrixPropagator::~LauKMatrixPropagator()
{
    // Destructor
    realProp_.Clear();
    negImagProp_.Clear();
    ScattKMatrix_.Clear();
    ReRhoMatrix_.Clear();
    ImRhoMatrix_.Clear();
    GammaMatrix_.Clear();
    ReTMatrix_.Clear();
    ImTMatrix_.Clear();
    IMatrix_.Clear();
    zeroMatrix_.Clear();
}

LauComplex LauKMatrixPropagator::getPropTerm( const Int_t channel ) const
{
    // Get the (i,j) = (index_, channel) term of the propagator
    // matrix. This allows us not to return the full propagator matrix.
    Double_t realTerm = this->getRealPropTerm( channel );
    Double_t imagTerm = this->getImagPropTerm( channel );

    LauComplex propTerm( realTerm, imagTerm );
    return propTerm;
}

Double_t LauKMatrixPropagator::getRealPropTerm( const Int_t channel ) const
{
    // Get the real part of the (i,j) = (index_, channel) term of the propagator
    // matrix. This allows us not to return the full propagator matrix.
    if ( parametersSet_ == kFALSE ) {
        return 0.0;
    }

    Double_t propTerm = realProp_[index_][channel];
    return propTerm;
}

Double_t LauKMatrixPropagator::getImagPropTerm( const Int_t channel ) const
{
    // Get the imaginary part of the (i,j) = (index_, channel) term of the propagator
    // matrix. This allows us not to return the full propagator matrix.
    if ( parametersSet_ == kFALSE ) {
        return 0.0;
    }

    Double_t imagTerm = -1.0 * negImagProp_[index_][channel];
    return imagTerm;
}

void LauKMatrixPropagator::updatePropagator( const Double_t s )
{
    // Calculate the K-matrix propagator for the given s value.
    // The K-matrix amplitude is given by
    // T_i = sum_{ij} (I - iK*rho)^-1 * P_j, where P is the production K-matrix.
    // i = index for the state (e.g. S-wave index = 0).
    // Here, we only find the (I - iK*rho)^-1 matrix part.

    // Check if we have almost the same s value as before. If so, don't re-calculate
    // the propagator nor any of the pole mass summation terms.
    if ( TMath::Abs( s - previousS_ ) < 1e-6 * s ) {
        //cout<<"Already got propagator for s = "<<s<<endl;
        return;
    }

    if ( parametersSet_ == kFALSE ) {
        //cerr<<"ERROR in LauKMatrixPropagator::updatePropagator. Parameters have not been set."<<endl;
        return;
    }

    // Calculate the denominator pole mass terms and Adler zero factor
    this->calcPoleDenomVect( s );
    this->updateAdlerZeroFactor( s );

    // Calculate the scattering K-matrix (real and symmetric)
    this->calcScattKMatrix( s );
    // Calculate the phase space density matrix, which is diagonal, but can be complex
    // if the quantity s is below various threshold values (analytic continuation).
    this->calcRhoMatrix( s );

    // Calculate the angular momentum barrier matrix, which is real and diagonal
    this->calcGammaMatrix( s );

    // Calculate K*rho*(gamma^2) (real and imaginary parts, since rho can be complex)
    TMatrixD GammaMatrixSq = ( GammaMatrix_ * GammaMatrix_ );
    TMatrixD K_realRhoGammaSq( ScattKMatrix_ );

    K_realRhoGammaSq *= ReRhoMatrix_;
    K_realRhoGammaSq *= GammaMatrixSq;
    TMatrixD K_imagRhoGammaSq( ScattKMatrix_ );

    K_imagRhoGammaSq *= ImRhoMatrix_;
    K_imagRhoGammaSq *= GammaMatrixSq;

    // A = I + K*Imag(rho)Gamma^2, B = -K*Real(Rho)Gamma^2
    // Calculate C and D matrices such that (A + iB)*(C + iD) = I,
    // ie. C + iD = (I - i K*rhoGamma^2)^-1, separated into real and imaginary parts.
    // realProp C = (A + B A^-1 B)^-1, imagProp D = -A^-1 B C
    TMatrixD A( IMatrix_ );
    A += K_imagRhoGammaSq;
    TMatrixD B( zeroMatrix_ );
    B -= K_realRhoGammaSq;

    TMatrixD invA( TMatrixD::kInverted, A );
    TMatrixD invA_B( invA );
    invA_B *= B;
    TMatrixD B_invA_B( B );
    B_invA_B *= invA_B;

    TMatrixD invC( A );
    invC += B_invA_B;
    // Set the real part of the propagator matrix ("C")
    realProp_ = TMatrixD( TMatrixD::kInverted, invC );

    // Set the (negative) imaginary part of the propagator matrix ("-D")
    TMatrixD BC( B );
    BC           *= realProp_;
    negImagProp_  = TMatrixD( invA );
    negImagProp_ *= BC;

    // Pre-multiply by the Gamma matrix:
    realProp_    = GammaMatrix_ * realProp_;
    negImagProp_ = GammaMatrix_ * negImagProp_;

    if ( verbose_ ) {
        std::cout << "In LauKMatrixPropagator::updatePropagator(s). D[1-iKrhoD^2]^-1: " << std::endl;
        TString realOutput( "Real part:" ), imagOutput( "Imag part:" );
        for ( int iChannel = 0; iChannel < nChannels_; iChannel++ ) {
            for ( int jChannel = 0; jChannel < nChannels_; jChannel++ ) {
                realOutput += Form( "\t%.6f", realProp_[iChannel][jChannel] );
                imagOutput += Form( "\t%.6f", -1 * negImagProp_[iChannel][jChannel] );
            }
            realOutput += "\n          ";
            imagOutput += "\n          ";
        }
        std::cout << realOutput << std::endl;
        std::cout << imagOutput << std::endl;
    }

    // Also calculate the production SVP term, since this uses Adler-zero parameters
    // defined in the parameter file.
    this->updateProdSVPTerm( s );

    // Finally, keep track of the value of s we just used.
    previousS_ = s;
}

void LauKMatrixPropagator::setParameters( const TString& inputFile )
{

    // Read an input file that specifies the values of the coupling constants g_i for
    // the given number of poles and their (bare) masses. Also provided are the f_{ab}
    // slow-varying constants. The input file should also provide the Adler zero
    // constants s_0, s_A and s_A0.

    parametersSet_ = kFALSE;

    cout << "Initialising K-matrix propagator " << name_ << " parameters from " << inputFile.Data()
         << endl;

    cout << "nChannels = " << nChannels_ << ", nPoles = " << nPoles_ << endl;

    // Initialise various matrices
    this->initialiseMatrices();
    std::vector<Int_t> a( nChannels_, 0 );

    // The format of the input file contains lines starting with a keyword followed by the
    // appropriate set of parameters. Keywords are case insensitive (treated as lower-case).
    // 1) Indices (nChannels) of N phase space channel types (defined in KMatrixChannels enum)
    // "Channels iChannel1 iChannel2 ... iChannelN"
    // 2) Definition of poles: bare mass (GeV), pole index (1 to NPoles), N channel couplings g_j
    // "Pole poleIndex mass g_1 g_2 ... g_N"
    // 3) Definition of scattering f_{ij} constants: scattering index (1 to N), channel values
    // "Scatt index f_{i1} f_{i2} ... f_{iN}", where i = index
    // 4) Orbital angular momentu for each channel. If not set here, defaults to 0
    // "AngularMomentum L[0] L[1] ... L[N]"
    // 5) Barrier factor parameter, which appears in the denominator and multiplies the term involving
    // the nominal radius. If not set here, defaults to 0 or values appropriate to non-zero angular
    // momenta as set in in (4) above.
    // "BarrierFactorParameter a[0] a[1] ... a[N]"
    // 6) Characteristic radius for each channel. If not set here, defaults to 3.0 GeV^{-1}
    // "Radii radChannel1 radChannel2 ... radChannelN"
    // 7) Various Adler zero and scattering constants; each line is "name value".
    // Possible names are mSq0, s0Scatt, s0Prod, sA, sA0
    //
    // By default, the scattering constants are symmetrised: f_{ji} = f_{ij}.
    // To not assume this use "ScattSymmetry 0" on a line

    LauTextFileParser readFile( inputFile );
    readFile.processFile();

    // Loop over the (non-commented) lines
    UInt_t nTotLines = readFile.getTotalNumLines();

    if ( nTotLines == 0 ) {
        std::cerr << "ERROR in LauKMatrixPropagator::setParameters : K-matrix parameter file not present - exiting."
                  << std::endl;

        gSystem->Exit( EXIT_FAILURE );
    }

    UInt_t iLine( 0 );

    for ( iLine = 1; iLine <= nTotLines; iLine++ ) {

        // Get the line of strings
        std::vector<std::string> theLine = readFile.getLine( iLine );

        // There should always be at least two strings: a keyword and at least 1 value
        if ( theLine.size() < 2 ) {
            continue;
        }

        TString keyword( theLine[0].c_str() );
        keyword.ToLower();    // Use lowercase

        if ( ! keyword.CompareTo( "channels" ) ) {

            // Channel indices for phase-space factors
            this->storeChannels( theLine );

        } else if ( ! keyword.CompareTo( "pole" ) ) {

            // Pole terms
            this->storePole( theLine );

        } else if ( ! keyword.CompareTo( "scatt" ) ) {

            // Scattering terms
            this->storeScattering( theLine );

        } else if ( ! keyword.CompareTo( "angularmomentum" ) ) {

            // Orbital angular momentum state for each channel & set default a values if called before storeBarrierFactorParameter
            this->storeOrbitalAngularMomenta( theLine, a );

        } else if ( ! keyword.CompareTo( "barrierfactorparameter" ) ) {

            // Value of parameter "a" in denominator of centrifugal barrier factor, gamma
            this->storeBarrierFactorParameter( theLine, a );

        } else if ( ! keyword.CompareTo( "radii" ) ) {

            // Values of characteristic radius
            this->storeRadii( theLine );

        } else {

            // Usually Adler-zero constants
            TString parString( theLine[1].c_str() );
            this->storeParameter( keyword, parString );
        }
    }

    sAConst_ = 0.5 * sA_.unblindValue() * LauConstants::mPiSq;

    // Symmetrise scattering parameters if enabled
    if ( scattSymmetry_ == kTRUE ) {

        for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {

            for ( Int_t jChannel = iChannel; jChannel < nChannels_; jChannel++ ) {

                LauParameter fPar                = fScattering_[iChannel][jChannel];
                fScattering_[jChannel][iChannel] = LauParameter( fPar );
            }
        }
    }

    // Now that radii and barrier-factor-denominator parameters have been set, cache the value of "a/(R*R)"
    for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {
        gamAInvRadSq_[iChannel] = a[iChannel] / ( radii_[iChannel] * radii_[iChannel] );
    }

    // All required parameters have been set
    parametersSet_ = kTRUE;

    cout << "Finished initialising K-matrix propagator " << name_ << endl;
}

void LauKMatrixPropagator::initialiseMatrices()
{

    // Identity and null matrices
    IMatrix_.Clear();
    IMatrix_.ResizeTo( nChannels_, nChannels_ );
    for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {
        IMatrix_[iChannel][iChannel] = 1.0;
    }

    zeroMatrix_.Clear();
    zeroMatrix_.ResizeTo( nChannels_, nChannels_ );

    // Real K matrix
    ScattKMatrix_.Clear();
    ScattKMatrix_.ResizeTo( nChannels_, nChannels_ );

    // Real and imaginary propagator matrices
    realProp_.Clear();
    negImagProp_.Clear();
    realProp_.ResizeTo( nChannels_, nChannels_ );
    negImagProp_.ResizeTo( nChannels_, nChannels_ );

    // Phase space matrices
    ReRhoMatrix_.Clear();
    ImRhoMatrix_.Clear();
    ReRhoMatrix_.ResizeTo( nChannels_, nChannels_ );
    ImRhoMatrix_.ResizeTo( nChannels_, nChannels_ );

    // Gamma matrices
    GammaMatrix_.Clear();
    GammaMatrix_.ResizeTo( nChannels_, nChannels_ );

    // Vector of orbital angular momenta for the channels (default is S-wave everywhere)
    L_.clear();
    L_.assign( nChannels_, 0 );

    // Characteristic radius (diagonal) vector (default to 3.0)
    radii_.clear();
    radii_.assign( nChannels_, 3.0 );

    // Vector to cache ratio a/R^2
    gamAInvRadSq_.clear();
    gamAInvRadSq_.resize( nChannels_ );

    // Square-root phase space matrices
    ReSqrtRhoMatrix_.Clear();
    ImSqrtRhoMatrix_.Clear();
    ReSqrtRhoMatrix_.ResizeTo( nChannels_, nChannels_ );
    ImSqrtRhoMatrix_.ResizeTo( nChannels_, nChannels_ );

    // T matrices
    ReTMatrix_.Clear();
    ImTMatrix_.Clear();
    ReTMatrix_.ResizeTo( nChannels_, nChannels_ );
    ImTMatrix_.ResizeTo( nChannels_, nChannels_ );

    // For the coupling and scattering constants, use LauParArrays instead of TMatrices
    // so that the quantities remain LauParameters instead of just doubles.
    // Each array is an stl vector of another stl vector of LauParameters:
    // std::vector< std::vector<LauParameter> >.
    // Set their sizes using the number of poles and channels defined in the constructor
    gCouplings_.clear();
    gCouplings_.resize( nPoles_ );

    for ( Int_t iPole = 0; iPole < nPoles_; iPole++ ) {
        gCouplings_[iPole].resize( nChannels_ );
    }

    fScattering_.clear();
    fScattering_.resize( nChannels_ );

    for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {
        fScattering_[iChannel].resize( nChannels_ );
    }

    // Clear other vectors
    phaseSpaceTypes_.clear();
    phaseSpaceTypes_.resize( nChannels_ );

    mSqPoles_.clear();
    mSqPoles_.resize( nPoles_ );

    haveCalled_storeBarrierFactorParameter = kFALSE;
}

void LauKMatrixPropagator::storeChannels( const std::vector<std::string>& theLine )
{

    // Get the list of channel indices to specify what phase space factors should be used
    // e.g. pipi, Kpi, eta eta', 4pi etc..

    // Check that the line has nChannels_+1 strings
    Int_t nTypes = static_cast<Int_t>( theLine.size() ) - 1;
    if ( nTypes != nChannels_ ) {
        cerr << "Error in LauKMatrixPropagator::storeChannels. The input file defines " << nTypes
             << " channels when " << nChannels_ << " are expected" << endl;
        return;
    }

    for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {

        Int_t phaseSpaceInt = std::atoi( theLine[iChannel + 1].c_str() );
        Bool_t checkChannel = this->checkPhaseSpaceType( phaseSpaceInt );

        if ( checkChannel == kTRUE ) {
            cout << "Adding phase space channel index " << phaseSpaceInt
                 << " to K-matrix propagator " << name_ << endl;
            phaseSpaceTypes_[iChannel] = static_cast<LauKMatrixPropagator::KMatrixChannels>(
                phaseSpaceInt );
        } else {
            cerr << "Phase space channel index " << phaseSpaceInt << " should be between 1 and "
                 << static_cast<int>( LauKMatrixPropagator::KMatrixChannels::TotChannels ) - 1
                 << endl;
        }
    }
}

void LauKMatrixPropagator::storePole( const std::vector<std::string>& theLine )
{

    // Store the pole mass and its coupling constants for each channel.
    // Each line will contain: Pole poleNumber poleMass poleCouplingsPerChannel

    // Check that the line has nChannels_ + 3 strings
    Int_t nWords  = static_cast<Int_t>( theLine.size() );
    Int_t nExpect = nChannels_ + 3;

    if ( nWords == nExpect ) {

        Int_t poleIndex = std::atoi( theLine[1].c_str() ) - 1;
        if ( poleIndex >= 0 && poleIndex < nPoles_ ) {

            Double_t poleMass   = std::atof( theLine[2].c_str() );
            Double_t poleMassSq = poleMass * poleMass;
            LauParameter mPoleParam( Form( "KM_%s_poleMassSq_%i", name_.Data(), poleIndex ),
                                     poleMassSq );
            mSqPoles_[poleIndex] = mPoleParam;

            cout << "Added bare pole mass " << poleMass << " GeV for pole number " << poleIndex + 1
                 << endl;

            for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {

                Double_t couplingConst = std::atof( theLine[iChannel + 3].c_str() );
                LauParameter couplingParam(
                    Form( "KM_%s_gCoupling_%i_%i", name_.Data(), poleIndex, iChannel ),
                    couplingConst );
                gCouplings_[poleIndex][iChannel] = couplingParam;

                cout << "Added coupling parameter g^{" << poleIndex + 1 << "}_" << iChannel + 1
                     << " = " << couplingConst << endl;
            }
        }

    } else {

        cerr << "Error in LauKMatrixPropagator::storePole. Expecting " << nExpect
             << " numbers for pole definition but found " << nWords << " values instead" << endl;
    }
}

void LauKMatrixPropagator::storeScattering( const std::vector<std::string>& theLine )
{

    // Store the scattering constants (along one of the channel rows).
    // Each line will contain: Scatt ScattIndex ScattConstantsPerChannel

    // Check that the line has nChannels_ + 2 strings
    Int_t nWords  = static_cast<Int_t>( theLine.size() );
    Int_t nExpect = nChannels_ + 2;

    if ( nWords == nExpect ) {

        Int_t scattIndex = std::atoi( theLine[1].c_str() ) - 1;
        if ( scattIndex >= 0 && scattIndex < nChannels_ ) {

            for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {

                Double_t scattConst = std::atof( theLine[iChannel + 2].c_str() );
                LauParameter scattParam(
                    Form( "KM_%s_fScatteringConst_%i_%i", name_.Data(), scattIndex, iChannel ),
                    scattConst );
                fScattering_[scattIndex][iChannel] = scattParam;

                cout << "Added scattering parameter f(" << scattIndex + 1 << "," << iChannel + 1
                     << ") = " << scattConst << endl;
            }
        }

    } else {

        cerr << "Error in LauKMatrixPropagator::storeScattering. Expecting " << nExpect
             << " numbers for scattering constants definition but found " << nWords
             << " values instead" << endl;
    }
}

void LauKMatrixPropagator::storeOrbitalAngularMomenta( const std::vector<std::string>& theLine,
                                                       std::vector<Int_t>& a )
{

    // Store the orbital angular momentum for each channel
    // Each line will contain: angularmomentum OrbitalAngularMomentumPerChannel

    // Check that the line has nChannels_ + 1 strings
    Int_t nWords  = static_cast<Int_t>( theLine.size() );
    Int_t nExpect = nChannels_ + 1;

    if ( nWords == nExpect ) {

        for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {

            Int_t angularMomentum = std::atoi( theLine[iChannel + 1].c_str() );
            L_[iChannel]          = angularMomentum;

            cout << "Defined K-matrix orbital angular momentum " << angularMomentum
                 << " for channel " << iChannel << endl;
        }

    } else {

        cerr << "Error in LauKMatrixPropagator::storeOrbitalAngularMomenta. Expecting " << nExpect
             << " numbers for orbital angular momenta definition but found " << nWords
             << " values instead" << endl;
    }

    if ( ! haveCalled_storeBarrierFactorParameter ) {
        // Set default value of spin-dependent centrifugal-barrier-factor parameter
        for ( Int_t iCh = 0; iCh < nChannels_; iCh++ ) {
            switch ( L_[iCh] ) {
                case 0 :
                    a[iCh] = 0;
                    break;
                case 1 :
                    a[iCh] = 1;
                    break;
                case 2 :
                    a[iCh] = 3;
                    break;
                default :
                    std::cerr << "ERROR in LauKMatrixPropagator constructor. Centrifugal barrier factor and angular-momentum terms of K-matrix are only defined for S-, P-, or D-wave."
                              << std::endl;
                    gSystem->Exit( EXIT_FAILURE );
            }
        }
    }
}

void LauKMatrixPropagator::storeRadii( const std::vector<std::string>& theLine )
{

    // Store the characteristic radii (measured in GeV^{-1})
    // Each line will contain: Radii RadiusConstantsPerChannel

    // Check that the line has nChannels_ + 1 strings
    Int_t nWords  = static_cast<Int_t>( theLine.size() );
    Int_t nExpect = nChannels_ + 1;

    if ( nWords == nExpect ) {

        for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {

            Double_t radiusConst = std::atof( theLine[iChannel + 1].c_str() );
            radii_[iChannel]     = radiusConst;

            cout << "Added K-matrix radius " << radiusConst << " for channel " << iChannel << endl;
        }

    } else {

        cerr << "Error in LauKMatrixPropagator::storeRadii. Expecting " << nExpect
             << " numbers for radii definition but found " << nWords << " values instead" << endl;
    }
}

void LauKMatrixPropagator::storeBarrierFactorParameter( const std::vector<std::string>& theLine,
                                                        std::vector<Int_t>& a )
{

    // Store the parameter of the barrier factor
    // Each line will contain: barrierfactorparameter ParameterValuePerchannel

    // Check that the line has nChannels_ + 1 strings
    Int_t nWords  = static_cast<Int_t>( theLine.size() );
    Int_t nExpect = nChannels_ + 1;

    if ( nWords == nExpect ) {

        for ( Int_t iChannel = 0; iChannel < nChannels_; iChannel++ ) {

            Double_t parameterValue = std::atof( theLine[iChannel + 1].c_str() );
            a[iChannel]             = parameterValue;

            cout << "Added K-matrix barrier factor parameter value " << parameterValue
                 << " for channel " << iChannel << endl;
        }

        // Set flag to stop storeOrbitalAngularMomenta overriding these a values
        haveCalled_storeBarrierFactorParameter = kTRUE;

    } else {

        cerr << "Error in LauKMatrixPropagator::storeBarrierFactorParameter. Expecting " << nExpect
             << " numbers for barrier factor parameter definition but found " << nWords
             << " values instead" << endl;
    }
}

void LauKMatrixPropagator::storeParameter( const TString& keyword, const TString& parString )
{

    if ( ! keyword.CompareTo( "msq0" ) ) {

        Double_t mSq0Value = std::atof( parString.Data() );
        cout << "Adler zero constant m0Sq = " << mSq0Value << endl;
        mSq0_ = LauParameter( Form( "KM_%s_mSq0", name_.Data() ), mSq0Value );

    } else if ( ! keyword.CompareTo( "s0scatt" ) ) {

        Double_t s0ScattValue = std::atof( parString.Data() );
        cout << "Adler zero constant s0Scatt = " << s0ScattValue << endl;
        s0Scatt_ = LauParameter( Form( "KM_%s_s0Scatt", name_.Data() ), s0ScattValue );

    } else if ( ! keyword.CompareTo( "s0prod" ) ) {

        Double_t s0ProdValue = std::atof( parString.Data() );
        cout << "Adler zero constant s0Prod = " << s0ProdValue << endl;
        s0Prod_ = LauParameter( Form( "KM_%s_s0Prod", name_.Data() ), s0ProdValue );

    } else if ( ! keyword.CompareTo( "sa0" ) ) {

        Double_t sA0Value = std::atof( parString.Data() );
        cout << "Adler zero constant sA0 = " << sA0Value << endl;
        sA0_ = LauParameter( Form( "KM_%s_sA0", name_.Data() ), sA0Value );

    } else if ( ! keyword.CompareTo( "sa" ) ) {

        Double_t sAValue = std::atof( parString.Data() );
        cout << "Adler zero constant sA = " << sAValue << endl;
        sA_ = LauParameter( Form( "KM_%s_sA", name_.Data() ), sAValue );

    } else if ( ! keyword.CompareTo( "scattsymmetry" ) ) {

        Int_t flag = std::atoi( parString.Data() );
        if ( flag == 0 ) {
            cout << "Turning off scattering parameter symmetry: f_ji = f_ij will not be assumed"
                 << endl;
            scattSymmetry_ = kFALSE;
        }
    }
}

void LauKMatrixPropagator::calcScattKMatrix( const Double_t s )
{

    // Calculate the scattering K-matrix for the given value of s.
    // We need to obtain the complete matrix (not just one row/column)
    // to get the correct inverted (I - i K rho) terms for the propagator.

    if ( verbose_ ) {
        cout << "Within calcScattKMatrix for s = " << s << endl;
    }

    // Initialise the K matrix to zero
    ScattKMatrix_.Zero();

    Int_t iChannel( 0 ), jChannel( 0 ), iPole( 0 );

    // The pole denominator 1/(m^2 - s) terms should already be calculated
    // by the calcPoleDenomVect() function. These same terms are also
    // used for calculating the production K-matrix elements.

    // Calculate the "slowly-varying part" (SVP), e.g. (1 GeV - s0)/(s - s0)
    this->updateScattSVPTerm( s );

    // Now loop over iChannel, jChannel to calculate Kij = Kji.
    for ( iChannel = 0; iChannel < nChannels_; iChannel++ ) {

        // Scattering matrix is real and symmetric. Start j loop from i.
        for ( jChannel = iChannel; jChannel < nChannels_; jChannel++ ) {

            Double_t Kij( 0.0 );

            // Calculate pole mass summation term
            for ( iPole = 0; iPole < nPoles_; iPole++ ) {

                Double_t g_i = this->getCouplingConstant( iPole, iChannel );
                Double_t g_j = this->getCouplingConstant( iPole, jChannel );

                Kij += poleDenomVect_[iPole] * g_i * g_j;
                if ( verbose_ ) {
                    cout << "1: Kij for i = " << iChannel << ", j = " << jChannel << " = " << Kij
                         << endl;
                }
            }

            Double_t fij  = this->getScatteringConstant( iChannel, jChannel );
            Kij          += fij * scattSVP_;

            Kij *= adlerZeroFactor_;
            if ( verbose_ ) {
                cout << "2: Kij for i = " << iChannel << ", j = " << jChannel << " = " << Kij << endl;
            }

            // Assign the TMatrix (i,j) element to the variable Kij and Kji (symmetry)
            ScattKMatrix_( iChannel, jChannel ) = Kij;
            ScattKMatrix_( jChannel, iChannel ) = Kij;

        }    // j loop

    }    // i loop
}

void LauKMatrixPropagator::calcPoleDenomVect( const Double_t s )
{
    // Calculate the 1/(m_pole^2 - s) terms for the scattering
    // and production K-matrix formulae.
    poleDenomVect_.clear();
    Int_t iPole( 0 );
    for ( iPole = 0; iPole < nPoles_; iPole++ ) {

        Double_t poleTerm = mSqPoles_[iPole].unblindValue() - s;
        Double_t invPoleTerm( 0.0 );
        if ( TMath::Abs( poleTerm ) > 1.0e-6 ) {
            invPoleTerm = 1.0 / poleTerm;
        }

        poleDenomVect_.push_back( invPoleTerm );
    }
}

Double_t LauKMatrixPropagator::getPoleDenomTerm( const Int_t poleIndex ) const
{

    if ( parametersSet_ == kFALSE ) {
        return 0.0;
    }

    Double_t poleDenom( 0.0 );
    poleDenom = poleDenomVect_[poleIndex];
    return poleDenom;
}

LauParameter& LauKMatrixPropagator::getPoleMassSqParameter( const Int_t poleIndex )
{
    if ( ( parametersSet_ == kFALSE ) || ( poleIndex < 0 || poleIndex >= nPoles_ ) ) {
        std::cerr << "ERROR from LauKMatrixPropagator::getPoleMassSqParameter(). Invalid pole."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    return mSqPoles_[poleIndex];
}

Double_t LauKMatrixPropagator::getCouplingConstant( const Int_t poleIndex,
                                                    const Int_t channelIndex ) const
{

    if ( parametersSet_ == kFALSE ) {
        return 0.0;
    }
    if ( poleIndex < 0 || poleIndex >= nPoles_ ) {
        return 0.0;
    }
    if ( channelIndex < 0 || channelIndex >= nChannels_ ) {
        return 0.0;
    }

    Double_t couplingConst = gCouplings_[poleIndex][channelIndex].unblindValue();
    return couplingConst;
}

LauParameter& LauKMatrixPropagator::getCouplingParameter( const Int_t poleIndex,
                                                          const Int_t channelIndex )
{

    if ( ( parametersSet_ == kFALSE ) || ( poleIndex < 0 || poleIndex >= nPoles_ ) ||
         ( channelIndex < 0 || channelIndex >= nChannels_ ) ) {
        std::cerr << "ERROR from LauKMatrixPropagator::getCouplingParameter(). Invalid coupling."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    //std::cout << "Minvalue + range for " << poleIndex << ", " << channelIndex << ": " << gCouplings_[poleIndex][channelIndex].minValue() << " => + " << gCouplings_[poleIndex][channelIndex].range() <<
    //             " and init value: " << gCouplings_[poleIndex][channelIndex].initValue() << std::endl;
    return gCouplings_[poleIndex][channelIndex];
}

Double_t LauKMatrixPropagator::getScatteringConstant( const Int_t channel1Index,
                                                      const Int_t channel2Index ) const
{

    if ( parametersSet_ == kFALSE ) {
        return 0.0;
    }
    if ( channel1Index < 0 || channel1Index >= nChannels_ ) {
        return 0.0;
    }
    if ( channel2Index < 0 || channel2Index >= nChannels_ ) {
        return 0.0;
    }

    Double_t scatteringConst = fScattering_[channel1Index][channel2Index].unblindValue();
    return scatteringConst;
}

LauParameter& LauKMatrixPropagator::getScatteringParameter( const Int_t channel1Index,
                                                            const Int_t channel2Index )
{

    if ( ( parametersSet_ == kFALSE ) || ( channel1Index < 0 || channel1Index >= nChannels_ ) ||
         ( channel2Index < 0 || channel2Index >= nChannels_ ) ) {
        std::cerr << "ERROR from LauKMatrixPropagator::getScatteringParameter(). Invalid chanel index."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    return fScattering_[channel1Index][channel2Index];
}

Double_t LauKMatrixPropagator::calcSVPTerm( const Double_t s, const Double_t s0 ) const
{

    if ( parametersSet_ == kFALSE ) {
        return 0.0;
    }

    // Calculate the "slowly-varying part" (SVP)
    Double_t result( 0.0 );
    Double_t deltaS = s - s0;
    if ( TMath::Abs( deltaS ) > 1.0e-6 ) {
        result = ( mSq0_.unblindValue() - s0 ) / deltaS;
    }

    return result;
}

void LauKMatrixPropagator::updateScattSVPTerm( const Double_t s )
{
    // Update the scattering "slowly-varying part" (SVP)
    Double_t s0Scatt = s0Scatt_.unblindValue();
    scattSVP_        = this->calcSVPTerm( s, s0Scatt );
}

void LauKMatrixPropagator::updateProdSVPTerm( const Double_t s )
{
    // Update the production "slowly-varying part" (SVP)
    Double_t s0Prod = s0Prod_.unblindValue();
    prodSVP_        = this->calcSVPTerm( s, s0Prod );
}

void LauKMatrixPropagator::updateAdlerZeroFactor( const Double_t s )
{

    // Calculate the multiplicative factor containing various Adler zero
    // constants.
    adlerZeroFactor_ = 0.0;

    Double_t sA0Val = sA0_.unblindValue();
    Double_t deltaS = s - sA0Val;
    if ( TMath::Abs( deltaS ) > 1e-6 ) {
        adlerZeroFactor_ = ( s - sAConst_ ) * ( 1.0 - sA0Val ) / deltaS;
    }
}

void LauKMatrixPropagator::calcGammaMatrix( const Double_t s )
{
    // Calculate the gamma angular momentum barrier matrix
    // for the given invariant mass squared quantity, s.

    // Initialise all entries to zero
    GammaMatrix_.Zero();

    Double_t gamma( 0.0 );

    for ( Int_t iChannel( 0 ); iChannel < nChannels_; ++iChannel ) {

        if ( L_[iChannel] != 0 ) {
            gamma = this->calcGamma( iChannel, s );
        } else {
            gamma = 1.0;    // S-wave
        }

        if ( verbose_ ) {
            cout << "GammaMatrix(" << iChannel << ", " << iChannel << ") = " << gamma << endl;
        }

        GammaMatrix_( iChannel, iChannel ) = gamma;
    }
}

Double_t LauKMatrixPropagator::calcGamma( const Int_t iCh, const Double_t s ) const
{
    // Calculate the barrier factor
    Double_t gamma( 0.0 );

    LauKMatrixPropagator::KMatrixChannels phaseSpaceIndex = phaseSpaceTypes_[iCh];

    LauComplex rho = getRho( s, phaseSpaceIndex );
    Double_t q     = 0.5 * sqrt( s ) * rho.abs();

    gamma = pow( q, L_[iCh] );
    if ( includeBWBarrierFactor_ ) {
        gamma /= pow( q * q + gamAInvRadSq_[iCh], L_[iCh] / 2. );
    }

    if ( verbose_ ) {
        std::cout << "In LauKMatrixPropagator::calcGamma(iCh=" << iCh << ", s=" << s << ", prop). ";
        std::cout << "|q(iCh=" << iCh << ")|: " << q << std::endl;
    }

    return gamma;
}

void LauKMatrixPropagator::calcRhoMatrix( const Double_t s )
{

    // Calculate the real and imaginary part of the phase space density
    // diagonal matrix for the given invariant mass squared quantity, s.
    // The matrix can be complex if s is below threshold (so that
    // the amplitude continues analytically).

    // Initialise all entries to zero
    ReRhoMatrix_.Zero();
    ImRhoMatrix_.Zero();

    for ( Int_t iChannel( 0 ); iChannel < nChannels_; ++iChannel ) {

        LauKMatrixPropagator::KMatrixChannels phaseSpaceIndex = phaseSpaceTypes_[iChannel];

        LauComplex rho = getRho( s, phaseSpaceIndex );

        if ( verbose_ ) {
            cout << "ReRhoMatrix(" << iChannel << ", " << iChannel << ") = " << rho.re() << endl;
            cout << "ImRhoMatrix(" << iChannel << ", " << iChannel << ") = " << rho.im() << endl;
        }

        ReRhoMatrix_( iChannel, iChannel ) = rho.re();
        ImRhoMatrix_( iChannel, iChannel ) = rho.im();
    }
}

LauComplex LauKMatrixPropagator::getRho( const Double_t s,
                                         const LauKMatrixPropagator::KMatrixChannels phaseSpaceIndex ) const
{
    LauComplex rho( 0.0, 0.0 );
    switch ( phaseSpaceIndex ) {
        case LauKMatrixPropagator::KMatrixChannels::PiPi :
            rho = this->calcPiPiRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::KK :
            rho = this->calcKKRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::FourPi :
            rho = this->calcFourPiRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::EtaEta :
            rho = this->calcEtaEtaRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::EtaEtaP :
            rho = this->calcEtaEtaPRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::KPi :
            rho = this->calcKPiRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::KEtaP :
            rho = this->calcKEtaPRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::KThreePi :
            rho = this->calcKThreePiRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::D0K :
            rho = this->calcD0KRho( s );
            break;
        case LauKMatrixPropagator::KMatrixChannels::Dstar0K :
            rho = this->calcDstar0KRho( s );
            break;
        default :
            std::cerr << "ERROR in LauKMatrixPropagator::getRho(...). Phase-space index not recognised for this channel"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
    }
    return rho;
}

LauComplex LauKMatrixPropagator::calcD0KRho( const Double_t s ) const
{
    // Calculate the D0K+ phase space factor
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm1 = ( -mD0KSumSq_ / s ) + 1.0;
    Double_t sqrtTerm2 = ( -mD0KDiffSq_ / s ) + 1.0;
    Double_t sqrtTerm  = sqrtTerm1 * sqrtTerm2;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcDstar0KRho( const Double_t s ) const
{
    // Calculate the Dstar0K+ phase space factor
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm1 = ( -mDstar0KSumSq_ / s ) + 1.0;
    Double_t sqrtTerm2 = ( -mDstar0KDiffSq_ / s ) + 1.0;
    Double_t sqrtTerm  = sqrtTerm1 * sqrtTerm2;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcPiPiRho( const Double_t s ) const
{
    // Calculate the pipi phase space factor
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm = ( -m2piSq_ / s ) + 1.0;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcKKRho( const Double_t s ) const
{
    // Calculate the KK phase space factor
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm = ( -m2KSq_ / s ) + 1.0;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcFourPiRho( const Double_t s ) const
{
    // Calculate the 4pi phase space factor. This uses a 6th-order polynomial
    // parameterisation that approximates the multi-body phase space double integral
    // defined in Eq 4 of the A&S paper hep-ph/0204328. This form agrees with the
    // BaBar model (another 6th order polynomial from s^4 down to 1/s^2), but avoids the
    // exponential increase at small values of s (~< 0.1) arising from 1/s and 1/s^2.
    // Eq 4 is evaluated for each value of s by assuming incremental steps of 1e-3 GeV^2
    // for s1 and s2, the invariant energy squared of each of the di-pion states,
    // with the integration limits of s1 = (2*mpi)^2 to (sqrt(s) - 2*mpi)^2 and
    // s2 = (2*mpi)^2 to (sqrt(s) - sqrt(s1))^2. The mass M of the rho is taken to be
    // 0.775 GeV and the energy-dependent width of the 4pi system
    // Gamma(s) = gamma_0*rho1^3(s), where rho1 = sqrt(1.0 - 4*mpiSq/s) and gamma_0 is
    // the "width" of the 4pi state at s = 1, which is taken to be 0.3 GeV
    // (~75% of the total width from PDG estimates of the f0(1370) -> 4pi state).
    // The normalisation term rho_0 is found by ensuring that the phase space integral
    // at s = 1 is equal to sqrt(1.0 - 16*mpiSq/s). Note that the exponent for this
    // factor in hep-ph/0204328 is wrong; it should be 0.5, i.e. sqrt, not n = 1 to 5.
    // Plotting the value of this double integral as a function of s can then be fitted
    // to a 6th-order polynomial (for s < 1), which is the result used below

    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    if ( s <= 1.0 ) {
        Double_t rhoTerm = ( ( 1.07885 * s + 0.13655 ) * s - 0.29744 ) * s - 0.20840;
        rhoTerm          = ( ( rhoTerm * s + 0.13851 ) * s - 0.01933 ) * s + 0.00051;
        // For some values of s (below 2*mpi), this term is a very small
        // negative number. Check for this and set the rho term to zero.
        if ( rhoTerm < 0.0 ) {
            rhoTerm = 0.0;
        }
        rho.setRealPart( rhoTerm );
    } else {
        rho.setRealPart( TMath::Sqrt( 1.0 - ( fourPiFactor1_ / s ) ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcEtaEtaRho( const Double_t s ) const
{
    // Calculate the eta-eta phase space factor
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm = ( -m2EtaSq_ / s ) + 1.0;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcEtaEtaPRho( const Double_t s ) const
{
    // Calculate the eta-eta' phase space factor. Note that the
    // mass difference term m_eta - m_eta' is not included,
    // since this corresponds to a "t or u-channel crossing",
    // which means that we cannot simply analytically continue
    // this part of the phase space factor below threshold, which
    // we can do for s-channel contributions. This is actually an
    // unsolved problem, e.g. see Guo et al 1409.8652, and
    // Danilkin et al 1409.7708. Anisovich and Sarantsev in
    // hep-ph/0204328 "solve" this issue by setting the mass
    // difference term to unity, which is what we do here...

    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm = ( -mEtaEtaPSumSq_ / s ) + 1.0;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcKPiRho( const Double_t s ) const
{
    // Calculate the K-pi phase space factor
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm1 = ( -mKpiSumSq_ / s ) + 1.0;
    Double_t sqrtTerm2 = ( -mKpiDiffSq_ / s ) + 1.0;
    Double_t sqrtTerm  = sqrtTerm1 * sqrtTerm2;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcKEtaPRho( const Double_t s ) const
{
    // Calculate the K-eta' phase space factor
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    Double_t sqrtTerm1 = ( -mKEtaPSumSq_ / s ) + 1.0;
    Double_t sqrtTerm2 = ( -mKEtaPDiffSq_ / s ) + 1.0;
    Double_t sqrtTerm  = sqrtTerm1 * sqrtTerm2;
    if ( sqrtTerm < 0.0 ) {
        rho.setImagPart( TMath::Sqrt( -sqrtTerm ) );
    } else {
        rho.setRealPart( TMath::Sqrt( sqrtTerm ) );
    }

    return rho;
}

LauComplex LauKMatrixPropagator::calcKThreePiRho( const Double_t s ) const
{
    // Calculate the Kpipipi + multimeson phase space factor.
    // Use the simplest definition in hep-ph/9705401 (Eq 14), which is the form
    // used for the rest of that paper (thankfully, the amplitude does not depend
    // significantly on the form used for the K3pi phase space factor).
    LauComplex rho( 0.0, 0.0 );
    if ( TMath::Abs( s ) < 1e-10 ) {
        return rho;
    }

    if ( s < 1.44 ) {

        Double_t powerTerm = ( -mK3piDiffSq_ / s ) + 1.0;
        if ( powerTerm < 0.0 ) {
            rho.setImagPart( k3piFactor_ * TMath::Power( -powerTerm, 2.5 ) );
        } else {
            rho.setRealPart( k3piFactor_ * TMath::Power( powerTerm, 2.5 ) );
        }

    } else {
        rho.setRealPart( 1.0 );
    }

    return rho;
}

Bool_t LauKMatrixPropagator::checkPhaseSpaceType( const Int_t phaseSpaceInt ) const
{
    Bool_t passed( kFALSE );

    if ( phaseSpaceInt >= 1 &&
         phaseSpaceInt < static_cast<Int_t>( LauKMatrixPropagator::KMatrixChannels::TotChannels ) ) {
        passed = kTRUE;
    }

    return passed;
}

LauComplex LauKMatrixPropagator::getTransitionAmp( const Double_t s, const Int_t channel )
{

    // Get the complex (unitary) transition amplitude T for the given channel
    LauComplex TAmp( 0.0, 0.0 );

    if ( channel <= 0 || channel > nChannels_ ) {
        return TAmp;
    }

    this->getTMatrix( s );

    TAmp.setRealPart( ReTMatrix_[index_][channel - 1] );
    TAmp.setImagPart( ImTMatrix_[index_][channel - 1] );

    return TAmp;
}

LauComplex LauKMatrixPropagator::getPhaseSpaceTerm( const Double_t s, const Int_t channel )
{

    // Get the complex (unitary) transition amplitude T for the given channel
    LauComplex rho( 0.0, 0.0 );

    if ( channel <= 0 || channel > nChannels_ ) {
        return rho;
    }

    // If s has changed from the previous value, recalculate rho
    if ( TMath::Abs( s - previousS_ ) > 1e-6 * s ) {
        this->calcRhoMatrix( s );
    }

    rho.setRealPart( ReRhoMatrix_[channel][channel - 1] );
    rho.setImagPart( ImRhoMatrix_[channel][channel - 1] );

    return rho;
}

void LauKMatrixPropagator::getTMatrix( const LauKinematics* kinematics )
{

    // Find the unitary T matrix, where T = [sqrt(rho)]^{*} T_hat sqrt(rho),
    // and T_hat = (I - i K rho)^-1 * K is the Lorentz-invariant T matrix,
    // which has phase-space factors included (rho). This function is not
    // needed to calculate the K-matrix amplitudes, but allows us
    // to check the variation of T as a function of s (kinematics)

    if ( ! kinematics ) {
        return;
    }

    // Get the invariant mass squared (s) from the kinematics object.
    // Use the resPairAmpInt to find which mass-squared combination to use.
    Double_t s( 0.0 );

    if ( resPairAmpInt_ == 1 ) {
        s = kinematics->getm23Sq();
    } else if ( resPairAmpInt_ == 2 ) {
        s = kinematics->getm13Sq();
    } else if ( resPairAmpInt_ == 3 ) {
        s = kinematics->getm12Sq();
    }

    this->getTMatrix( s );
}

void LauKMatrixPropagator::getTMatrix( const Double_t s )
{

    // Find the unitary transition T matrix, where
    // T = [sqrt(rho)]^{*} T_hat sqrt(rho), and
    // T_hat = (I - i K rho)^-1 * K is the Lorentz-invariant T matrix,
    // which has phase-space factors included (rho). Note that the first
    // sqrt of the rho matrix is complex conjugated.

    // This function is not needed to calculate the K-matrix amplitudes, but
    // allows us to check the variation of T as a function of s (kinematics)

    // Initialse the real and imaginary parts of the T matrix to zero
    ReTMatrix_.Zero();
    ImTMatrix_.Zero();

    if ( parametersSet_ == kFALSE ) {
        return;
    }

    // Update K, rho and the propagator (I - i K rho)^-1
    this->updatePropagator( s );

    // Find the real and imaginary T_hat matrices
    TMatrixD THatReal = realProp_ * ScattKMatrix_;
    TMatrixD THatImag( zeroMatrix_ );
    THatImag -= negImagProp_ * ScattKMatrix_;

    // Find the square-root of the phase space matrix
    this->getSqrtRhoMatrix();

    // Let sqrt(rho) = A + iB and T_hat = C + iD
    // => T = A(CA-DB) + B(DA+CB) + i[A(DA+CB) + B(DB-CA)]
    TMatrixD CA( THatReal );
    CA *= ReSqrtRhoMatrix_;
    TMatrixD DA( THatImag );
    DA *= ReSqrtRhoMatrix_;
    TMatrixD CB( THatReal );
    CB *= ImSqrtRhoMatrix_;
    TMatrixD DB( THatImag );
    DB *= ImSqrtRhoMatrix_;

    TMatrixD CAmDB( CA );
    CAmDB -= DB;
    TMatrixD DApCB( DA );
    DApCB += CB;
    TMatrixD DBmCA( DB );
    DBmCA -= CA;

    // Find the real and imaginary parts of the transition matrix T
    ReTMatrix_ = ReSqrtRhoMatrix_ * CAmDB + ImSqrtRhoMatrix_ * DApCB;
    ImTMatrix_ = ReSqrtRhoMatrix_ * DApCB + ImSqrtRhoMatrix_ * DBmCA;
}

void LauKMatrixPropagator::getSqrtRhoMatrix()
{

    // Find the square root of the (current) phase space matrix so that
    // we can find T = [sqrt(rho)}^{*} T_hat sqrt(rho), where T_hat is the
    // Lorentz-invariant T matrix = (I - i K rho)^-1 * K; note that the first
    // sqrt of rho matrix is complex conjugated

    // If rho = rho_i + i rho_r = a + i b, then sqrt(rho) = c + i d, where
    // c = sqrt(0.5*(r+a)) and d = sqrt(0.5(r-a)), where r = sqrt(a^2 + b^2).
    // Since rho is diagonal, then the square root of rho will also be diagonal,
    // with its real and imaginary matrix elements equal to c and d, respectively

    // Initialise the real and imaginary parts of the square root of
    // the rho matrix to zero
    ReSqrtRhoMatrix_.Zero();
    ImSqrtRhoMatrix_.Zero();

    for ( Int_t iChannel( 0 ); iChannel < nChannels_; ++iChannel ) {

        Double_t realRho = ReRhoMatrix_[iChannel][iChannel];
        Double_t imagRho = ImRhoMatrix_[iChannel][iChannel];

        Double_t rhoMag  = sqrt( realRho * realRho + imagRho * imagRho );
        Double_t rhoSum  = rhoMag + realRho;
        Double_t rhoDiff = rhoMag - realRho;

        Double_t reSqrtRho( 0.0 ), imSqrtRho( 0.0 );
        if ( rhoSum > 0.0 ) {
            reSqrtRho = sqrt( 0.5 * rhoSum );
        }
        if ( rhoDiff > 0.0 ) {
            imSqrtRho = sqrt( 0.5 * rhoDiff );
        }

        ReSqrtRhoMatrix_[iChannel][iChannel] = reSqrtRho;
        ImSqrtRhoMatrix_[iChannel][iChannel] = imSqrtRho;
    }
}

LauComplex LauKMatrixPropagator::getTHat( const Double_t s, const Int_t channel )
{

    LauComplex THat( 0.0, 0.0 );

    if ( channel <= 0 || channel > nChannels_ ) {
        return THat;
    }

    this->updatePropagator( s );

    // Find the real and imaginary T_hat matrices
    TMatrixD THatReal = realProp_ * ScattKMatrix_;
    TMatrixD THatImag( zeroMatrix_ );
    THatImag -= negImagProp_ * ScattKMatrix_;

    // Return the specific THat component
    THat.setRealPart( THatReal[index_][channel - 1] );
    THat.setImagPart( THatImag[index_][channel - 1] );

    return THat;
}
