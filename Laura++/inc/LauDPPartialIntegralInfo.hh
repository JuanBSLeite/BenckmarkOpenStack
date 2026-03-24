
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

/*! \file LauDPPartialIntegralInfo.hh
    \brief File containing declaration of LauDPPartialIntegralInfo class.
*/

/*! \class LauDPPartialIntegralInfo
    \brief Class for defining (a section of) the Dalitz plot integration binning scheme

    Defines the range and bin size of the integration grid.
    Stores the weights and Jacobian for each grid point.
    Also stores the amplitude values for each model component at each grid point.
*/

#ifndef LAU_DPPARTIALINTEGRAL_INFO
#define LAU_DPPARTIALINTEGRAL_INFO

#include "LauComplex.hh"

#include "TString.h"

#include <iosfwd>

class LauKinematics;

class LauDPPartialIntegralInfo {

  public:
    //! Constructor
    /*!
        \param [in] minm13 the minimum of the m13 range
        \param [in] maxm13 the maximum of the m13 range
        \param [in] minm23 the minimum of the m23 range
        \param [in] maxm23 the maximum of the m23 range
        \param [in] m13BinWidth the m13 bin width
        \param [in] m23BinWidth the m23 bin width
        \param [in] precision the precision required for the Gauss-Legendre weights
        \param [in] nAmp the number of coherent amplitude components
        \param [in] nIncohAmp the number of incoherent amplitude components
        \param [in] squareDP whether or not to use the square DP for the integration - if so, m13 is actually mPrime and m23 is actually thetaPrime
        \param [in] kinematics the kinematics object to use to calculate the Jacobians (only relevant if squareDP is true)
    */
    LauDPPartialIntegralInfo( const Double_t minm13,
                              const Double_t maxm13,
                              const Double_t minm23,
                              const Double_t maxm23,
                              const Double_t m13BinWidth,
                              const Double_t m23BinWidth,
                              const Double_t precision,
                              const UInt_t nAmp,
                              const UInt_t nIncohAmp,
                              const Bool_t squareDP           = kFALSE,
                              const LauKinematics* kinematics = 0 );

    //! Destructor
    virtual ~LauDPPartialIntegralInfo();

    //! Retrieve the minm13 of DP
    /*!
        \return the the minm13 of DP
    */
    inline Double_t getMinm13() const { return minm13_; }

    //! Retrieve the maxm13 of DP
    /*!
        \return the the maxm13 of DP
    */
    inline Double_t getMaxm13() const { return maxm13_; }

    //! Retrieve the minm23 of DP
    /*!
        \return the the minm23 of DP
    */
    inline Double_t getMinm23() const { return minm23_; }

    //! Retrieve the maxm23 of DP
    /*!
        \return the the maxm23 of DP
    */
    inline Double_t getMaxm23() const { return maxm23_; }

    //! Retrieve the m13BinWidth of DP
    /*!
        \return the the m13BinWidth of DP
    */
    inline Double_t getM13BinWidth() const { return m13BinWidth_; }

    //! Retrieve the m23BinWidth of DP
    /*!
        \return the the m23BinWidth of DP
    */
    inline Double_t getM23BinWidth() const { return m23BinWidth_; }

    //! Retrieve the number of bins in m13
    /*!
        \return the number of bins in m13
    */
    inline UInt_t getnm13Points() const { return nm13Points_; }

    //! Retrieve the number of bins in m23
    /*!
        \return the number of bins in m23
    */
    inline UInt_t getnm23Points() const { return nm23Points_; }

    //! Retrieve the square DP flag
    /*!
        \return whether or not the integration is performed in the square DP
    */
    inline Bool_t getSquareDP() const { return squareDP_; }

    //! Retrieve the weight for the given grid point
    /*!
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
        \return the value of the weight
    */
    inline Double_t getWeight( const UInt_t m13Point, const UInt_t m23Point ) const
    {
        return weights_[m13Point][m23Point];
    }

    //! Retrieve the m13 value at the given grid point
    /*!
        \param [in] m13Point the grid index in m13
        \return the m13 value
    */
    inline Double_t getM13Value( const UInt_t m13Point ) const { return m13Points_[m13Point]; }

    //! Retrieve the m23 value at the given grid point
    /*!
        \param [in] m23Point the grid index in m23
        \return the m23 value
    */
    inline Double_t getM23Value( const UInt_t m23Point ) const { return m23Points_[m23Point]; }

    //! Retrieve the efficiency for the given grid point
    /*!
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
        \return the efficiency value
    */
    inline Double_t getEfficiency( const UInt_t m13Point, const UInt_t m23Point ) const
    {
        return efficiencies_[m13Point][m23Point];
    }

    //! Store the efficiency for the given grid point
    /*!
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
        \param [in] efficiency the new efficiency value
    */
    inline void storeEfficiency( const UInt_t m13Point,
                                 const UInt_t m23Point,
                                 const Double_t efficiency )
    {
        efficiencies_[m13Point][m23Point] = efficiency;
    }

    //! Retrieve the amplitude for the given grid point and amplitude index
    /*!
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
        \param [in] iAmp the amplitude index
        \return the amplitude value
    */
    inline const LauComplex& getAmplitude( const UInt_t m13Point,
                                           const UInt_t m23Point,
                                           const UInt_t iAmp ) const
    {
        return amplitudes_[m13Point][m23Point][iAmp];
    }

    //! Store the amplitude for the given grid point and amplitude index
    /*!
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
        \param [in] iAmp the amplitude index
        \param [in] amplitude the new amplitude value
    */
    inline void storeAmplitude( const UInt_t m13Point,
                                const UInt_t m23Point,
                                const UInt_t iAmp,
                                const LauComplex& amplitude )
    {
        amplitudes_[m13Point][m23Point][iAmp] = amplitude;
    }

    //! Retrieve the intensity for the given grid point and intensity index
    /*!
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
        \param [in] iAmp the intensity index
        \return the intensity value
    */
    inline Double_t getIntensity( const UInt_t m13Point, const UInt_t m23Point, const UInt_t iAmp ) const
    {
        return incohIntensities_[m13Point][m23Point][iAmp];
    }

    //! Store the intensity for the given grid point and intensity index
    /*!
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
        \param [in] iAmp the intensity index
        \param [in] intensity the new intensity value
    */
    inline void storeIntensity( const UInt_t m13Point,
                                const UInt_t m23Point,
                                const UInt_t iAmp,
                                const Double_t intensity )
    {
        incohIntensities_[m13Point][m23Point][iAmp] = intensity;
    }

  private:
    //! Copy constructor (not implemented)
    LauDPPartialIntegralInfo( const LauDPPartialIntegralInfo& other );

    //! Copy assignment operator (not implemented)
    LauDPPartialIntegralInfo& operator=( const LauDPPartialIntegralInfo& other );

    //! The minimum of the m13 range
    const Double_t minm13_;

    //! The maximum of the m13 range
    const Double_t maxm13_;

    //! The minimum of the m23 range
    const Double_t minm23_;

    //! The maximum of the m23 range
    const Double_t maxm23_;

    //! The bin width for m13
    const Double_t m13BinWidth_;

    //! The bin width for m23
    const Double_t m23BinWidth_;

    //! The number of bins in m13
    const UInt_t nm13Points_;

    //! The number of bins in m23
    const UInt_t nm23Points_;

    //! The number of amplitude components
    const UInt_t nAmp_;

    //! The number of amplitude components
    const UInt_t nIncohAmp_;

    //! Flag whether or not we're using the square DP for the integration
    const Bool_t squareDP_;

    //! The m13 positions of the grid points
    std::vector<Double_t> m13Points_;

    //! The m23 positions of the grid points
    std::vector<Double_t> m23Points_;

    //! The Gauss-Legendre weights of the m13 grid points
    std::vector<Double_t> m13Weights_;

    //! The Gauss-Legendre weights of the m23 grid points
    std::vector<Double_t> m23Weights_;

    //! The combined weights at each 2D grid point
    std::vector<std::vector<Double_t>> weights_;

    //! The efficiency at each 2D grid point
    std::vector<std::vector<Double_t>> efficiencies_;

    //! The amplitude values at each 2D grid point
    std::vector<std::vector<std::vector<LauComplex>>> amplitudes_;

    //! The incoherent intensity values at each 2D grid point
    std::vector<std::vector<std::vector<Double_t>>> incohIntensities_;

    ClassDef( LauDPPartialIntegralInfo, 0 )
};

//! output operator formatting
std::ostream& operator<<( std::ostream& stream, const LauDPPartialIntegralInfo& infoRecord );

#endif
