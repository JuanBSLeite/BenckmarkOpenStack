
/*
Copyright 2006 University of Warwick

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

/*! \file LauScfMap.hh
    \brief File containing declaration of LauScfMap class.
*/

/*! \class LauScfMap
    \brief Class for representing the 4D smearing matrix for mis-reconstructed signal (self cross feed)

    The smearing matrix is formed from a 2D array of 2D histograms
*/

#ifndef LAU_SCF_MAP
#define LAU_SCF_MAP

#include "Rtypes.h"

#include <map>
#include <vector>

class TH2;

class LauScfMap {

  public:
    //! Constructor
    LauScfMap();

    //! Destructor
    virtual ~LauScfMap();

    //! Construct the smearing matrix from the supplied histograms
    /*!
        The histograms should give the migration probabilities from a given true bin to a given reco bin.

        These can be constructed as follows:\n
        Each histogram contains only the events that were generated in a given true bin.\n
        The events are plotted at their reconstructed co-ordinates.\n
        The histogram is then normalised.

        All histograms should be identically binned.\n
        There should be as many histograms supplied as there are bins in each.\n
        The order of the histograms in the vector should match the ROOT "global bin number" ordering.

        \param [in] histos vector of two-dimensional histograms
    */
    void setHistos( const std::vector<const TH2*>& histos );

    //! Get the number of bins
    /*!
        \return the number of bins
    */
    Int_t getNumberOfBins() const { return numberOfBins_; };

    //! Create lists of the co-ordinates of the centres of true bins
    /*!
        \param [out] xCoords the vector of x co-ordinates of the centre of each true bin
        \param [out] yCoords the vector of y co-ordinates of the centre of each true bin
    */
    void listBinCentres( std::vector<Double_t>& xCoords, std::vector<Double_t>& yCoords ) const;

    //! Find the global bin number for the given co-ordinates
    /*!
        \param [in] xCoord the vector of x co-ordinate
        \param [in] yCoord the vector of y co-ordinate
        \return the global bin number
    */
    Int_t binNumber( Double_t xCoord, Double_t yCoord ) const;

    //! Find which true bins contribute to the given reco bin
    /*!
        \param [in] recoBin the reco global bin number
        \return the list of true bins that contribute (the numbers returned are not global bin numbers but the index in the LauScfMap::histos_ storage)
    */
    const std::vector<Int_t>* trueBins( Int_t recoBin ) const;

    //! Probability of a true event in the given true bin migrating to the reco bin
    /*!
        \param [in] recoBin the reco global bin number
        \param [in] trueBin the true bin index in the LauScfMap::histos_ storage (as provided by LauScfMap::trueBins)
        \return the migration probability
    */
    Double_t prob( Int_t recoBin, Int_t trueBin ) const;

    //! Retrieve the migration histogram for trueBin
    /*!
        \param [in] trueBin the true global bin number
        \return the migration histogram
    */
    TH2* trueHist( Int_t trueBin );

  private:
    //! Copy constructor (not implemented)
    LauScfMap( const LauScfMap& rhs );

    //! Copy assignment operator (not implemented)
    LauScfMap& operator=( const LauScfMap& rhs );

    //! The vector of two-dimensional histograms
    std::vector<TH2*> histos_;

    //! Map that links each reco global bin number and the vector of true bin indices into LauScfMap::histos_
    std::map<Int_t, std::vector<Int_t>> contribs_;

    //! Number of bins in the x-axis
    Int_t nBinsX_;
    //! Number of bins in the y-axis
    Int_t nBinsY_;
    //! Number of bins in a 2D histogram (nBinsX_*nBinsY_)
    Int_t numberOfBins_;

    ClassDef( LauScfMap, 0 )
};

#endif
