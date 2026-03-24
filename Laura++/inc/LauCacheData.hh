
/*
Copyright 2005 University of Warwick

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

/*! \file LauCacheData.hh
    \brief File containing declaration of LauCacheData class.
*/

/*! \class LauCacheData
    \brief Class to contain cached data relating to an event.

    Contains information on the dynamics of the event and some additional information such as the efficiency and the self cross feed fraction.
    Used by the DP dynamics classes to store the cached event information.
*/

#ifndef LAU_CACHE_DATA
#define LAU_CACHE_DATA

#include <vector>

#include <Rtypes.h>

class LauCacheData {

  public:
    //! Constructor
    LauCacheData();

    //! Destructor
    virtual ~LauCacheData();

    //! Copy constructor
    /*!
        \param [in] rhs the object to copy
    */
    LauCacheData( const LauCacheData& rhs );

    //! Copy assignment operator
    /*!
        \param [in] rhs the object to clone
    */
    LauCacheData& operator=( const LauCacheData& rhs );

    //! Set the invariant mass squared of the first and third daugthers
    /*!
        \param [in] m13Sq the invariant mass squared of the first and third daugthers
    */
    inline void storem13Sq( Double_t m13Sq ) { m13Sq_ = m13Sq; }

    //! Set the invariant mass squared of the second and third daugthers
    /*!
        \param [in] m23Sq the invariant mass squared of the second and third daugthers
    */
    inline void storem23Sq( Double_t m23Sq ) { m23Sq_ = m23Sq; }

    //! Set the square Dalitz plot coordinate, m'
    /*!
        \param [in] mPrime the square Dalitz plot coordinate, m'
    */
    inline void storemPrime( Double_t mPrime ) { mPrime_ = mPrime; }

    //! Set the square Dalitz plot coordinate, theta'
    /*!
        \param [in] thPrime the square Dalitz plot coordinate, theta'
    */
    inline void storethPrime( Double_t thPrime ) { thPrime_ = thPrime; }

    //! Set the tagging category
    /*!
        \param [in] tagCat the tagging category
    */
    inline void storeTagCat( Int_t tagCat ) { tagCat_ = tagCat; }

    //! Set the efficiency
    /*!
        \param [in] eff the efficiency
    */
    inline void storeEff( Double_t eff ) { eff_ = eff; }

    //! Set the fraction of poorly constructed events (the self cross feed fraction)
    /*!
        \param [in] scfFraction the fraction of poorly constructed events
    */
    inline void storeScfFraction( Double_t scfFraction ) { scfFraction_ = scfFraction; }

    //! Set the Jacobian for the transformation into square Dalitz coordinates
    /*!
        \param [in] jacobian the Jacobian
    */
    inline void storeJacobian( Double_t jacobian ) { jacobian_ = jacobian; }

    //! Set the real parts of the amplitudes
    /*!
        \param [in] realAmp the real parts of the amplitudes
    */
    inline void storeRealAmp( const std::vector<Double_t>& realAmp ) { realAmp_ = realAmp; }

    //! Set the imaginary parts of the amplitudes
    /*!
        \param [in] imagAmp the imaginary parts of the amplitudes
    */
    inline void storeImagAmp( const std::vector<Double_t>& imagAmp ) { imagAmp_ = imagAmp; }

    //! Set the incoherent intensities
    /*!
        \param [in] incohIntensities the intensities of the incoherent contributions
    */
    inline void storeIncohIntensities( const std::vector<Double_t>& incohIntensities )
    {
        incohIntensities_ = incohIntensities;
    }

    //! Retrieve the invariant mass squared of the first and third daugthers
    /*!
        \return m13 squared
    */
    inline Double_t retrievem13Sq() const { return m13Sq_; }

    //! Retrieve the invariant mass squared of the second and third daugthers
    /*!
        \return m23 squared
    */
    inline Double_t retrievem23Sq() const { return m23Sq_; }

    //! Retrieve the square Dalitz plot coordinate, m'
    /*!
        \return m'
    */
    inline Double_t retrievemPrime() const { return mPrime_; }

    //! Retrieve the square Dalitz plot coordinate, theta'
    /*!
        \return theta'
    */
    inline Double_t retrievethPrime() const { return thPrime_; }

    //! Retrieve the tagging category
    /*!
        \return the tagging category
    */
    inline Int_t retrieveTagCat() const { return tagCat_; }

    //! Retrieve the efficiency
    /*!
        \return the efficiency
    */
    inline Double_t retrieveEff() const { return eff_; }

    //! Retrieve the fraction of poorly reconstructed events (the so-called self cross feed fraction)
    /*!
        \return the fraction of poorly reconstructed events
    */
    inline Double_t retrieveScfFraction() const { return scfFraction_; }

    //! Retrieve the Jacobian for the transformation into square-Dalitz-plot coordinates
    /*!
        \return the Jacobian
    */
    inline Double_t retrieveJacobian() const { return jacobian_; }

    //! Retrieve the real parts of the amplitudes
    /*!
        \return the real parts of the amplitudes
    */
    inline const std::vector<Double_t>& retrieveRealAmp() const { return realAmp_; }

    //! Retrieve the imaginary parts of the amplitudes
    /*!
        \return the imaginary parts of the amplitudes
    */
    inline const std::vector<Double_t>& retrieveImagAmp() const { return imagAmp_; }

    //! Retrieve the incoherent intensities
    /*!
        \return the intensities of the incoherent contributions
    */
    inline const std::vector<Double_t>& retrieveIncohIntensities() const
    {
        return incohIntensities_;
    }

    //! Retrieve the real parts of the amplitudes
    /*!
        \return the real parts of the amplitudes
    */
    inline std::vector<Double_t>& retrieveRealAmp() { return realAmp_; }

    //! Retrieve the imaginary parts of the amplitudes
    /*!
        \return the imaginary parts of the amplitudes
    */
    inline std::vector<Double_t>& retrieveImagAmp() { return imagAmp_; }

    //! Retrieve the incoherent intensities
    /*!
        \return the intensities of the incoherent contributions
    */
    inline std::vector<Double_t>& retrieveIncohIntensities() { return incohIntensities_; }

  protected:

  private:
    //! The invariant mass squared of the first and third daugthers
    Double_t m13Sq_;

    //! The invariant mass squared of the second and third daugthers
    Double_t m23Sq_;

    //! The square Dalitz plot coordinate, m'
    Double_t mPrime_;

    //! The square Dalitz plot coordinate, theta'
    Double_t thPrime_;

    //! The tagging category
    Int_t tagCat_;

    //! The efficiency
    Double_t eff_;

    //! The fraction of poorly constructed events (the self cross feed fraction)
    Double_t scfFraction_;

    //! The Jacobian for the transformation into square Dalitz coordinates
    Double_t jacobian_;

    //! The real parts of the amplitudes
    std::vector<Double_t> realAmp_;

    //! The imaginary parts of the amplitudes
    std::vector<Double_t> imagAmp_;

    //! The intensities of the incoherent contributions
    std::vector<Double_t> incohIntensities_;

    ClassDef( LauCacheData, 0 )    // Cached Data Class
};

#endif
