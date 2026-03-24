
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

/*! \file LauComplex.hh
    \brief File containing declaration of LauComplex class.
*/

/*! \class LauComplex
    \brief Class for defining a complex number

    Class for complex number manipulation.
    In the function descriptions, the form (a,b) is used to represent a complex
    number. This is equivalent to the mathematical expression a + ib.
*/

/*****************************************************************************
 * Class based on RooFit/RooComplex.                                         *
 * Original copyright given below.                                           *
 *****************************************************************************
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef LAU_COMPLEX
#define LAU_COMPLEX

#include "Rtypes.h"
#include "TMath.h"

#include <iosfwd>

class LauComplex {

  public:
    //! Default Constructor
    inline LauComplex() :
        re_( 0.0 ),
        im_( 0.0 )
    {
    }

    //! Constructor
    /*!
        \param [in] a the value corresponding to the real part of the complex number
        \param [in] b the value corresponding to the imaginary part of the complex number
    */
    inline LauComplex( Double_t a, Double_t b ) :
        re_( a ),
        im_( b )
    {
    }

    //! Destructor
    virtual ~LauComplex() {}

    //! Copy constructor
    /*!
        \param [in] other the complex number to be copied
    */
    inline LauComplex( const LauComplex& other ) :
        re_( other.re_ ),
        im_( other.im_ )
    {
    }

    //! Copy assignment operator
    /*!
        \param [in] other the complex number to be copied
        \return the assigned complex number
    */
    inline LauComplex& operator=( const LauComplex& other )
    {
        if ( &other != this ) {
            re_ = other.re_;
            im_ = other.im_;
        }
        return *this;
    }

    //! Unary minus operator
    /*!
        \return the negated complex number
    */
    inline LauComplex operator-() const { return LauComplex( -re_, -im_ ); }

    //! Addition operator
    /*!
        \param [in] other the object to added to this one
        \return the sum of the two complex numbers
    */
    inline LauComplex operator+( const LauComplex& other ) const
    {
        LauComplex tmpCmplx( *this );
        tmpCmplx += other;
        return tmpCmplx;
    }

    //! Subtraction operator
    /*!
        \param [in] other the object to be subtracted from this one
        \return the difference of the two complex numbers
    */
    inline LauComplex operator-( const LauComplex& other ) const
    {
        LauComplex tmpCmplx( *this );
        tmpCmplx -= other;
        return tmpCmplx;
    }

    //! Multiplication operator
    /*!
        \param [in] other the object this one is to be multiplied by
        \return the product of the two complex numbers
    */
    inline LauComplex operator*( const LauComplex& other ) const
    {
        LauComplex tmpCmplx( *this );
        tmpCmplx *= other;
        return tmpCmplx;
    }

    //! Division operator
    /*!
        \param [in] other the object this one is to be divided by
        \return the ratio of the two complex numbers
    */
    inline LauComplex operator/( const LauComplex& other ) const
    {
        LauComplex tmpCmplx( *this );
        tmpCmplx /= other;
        return tmpCmplx;
    }

    //! Addition assignment operator
    /*!
        \param [in] other the object to be added to this one
        \return the result of the addition
    */
    inline LauComplex operator+=( const LauComplex& other )
    {
        this->re_ += other.re_;
        this->im_ += other.im_;
        return ( *this );
    }

    //! Subtraction assignment operator
    /*!
        \param [in] other the object to be subtracted from this one
        \return the result of the subtraction
    */
    inline LauComplex operator-=( const LauComplex& other )
    {
        this->re_ -= other.re_;
        this->im_ -= other.im_;
        return ( *this );
    }

    //! Multiplication assignment operator
    /*!
        \param [in] other the object this one is to be multiplied by
        \return the result of the multiplication
    */
    inline LauComplex operator*=( const LauComplex& other )
    {
        Double_t realPart = this->re_ * other.re_ - this->im_ * other.im_;
        Double_t imagPart = this->re_ * other.im_ + this->im_ * other.re_;
        this->setRealImagPart( realPart, imagPart );
        return ( *this );
    }

    //! Division assignment operator
    /*!
        \param [in] other the object this one is to be divided by
        \return the results of the division
    */
    inline LauComplex operator/=( const LauComplex& other )
    {
        Double_t x( other.abs2() );
        Double_t realPart = ( this->re_ * other.re_ + this->im_ * other.im_ ) / x;
        Double_t imagPart = ( this->im_ * other.re_ - this->re_ * other.im_ ) / x;
        this->setRealImagPart( realPart, imagPart );
        return ( *this );
    }

    //! Boolean comparison operator
    /*!
        \param [in] other the object to compared with this one
        \return true/false for the comparison
    */
    inline Bool_t operator==( const LauComplex& other ) const
    {
        return ( re_ == other.re_ && im_ == other.im_ );
    }

    //! Get the real part
    /*!
        \return the real part of the complex number
    */
    inline Double_t re() const { return re_; }

    //! Get the imaginary part
    /*!
        \return the imaginary part of the complex number
    */
    inline Double_t im() const { return im_; }

    //! Obtain the absolute value of the complex number
    /*!
        \return the absolute value (magnitude)
    */
    inline Double_t abs() const { return TMath::Sqrt( this->abs2() ); }

    //! Obtain the square of the absolute value of the complex number
    /*!
        \return the square of the absolute value (magnitude^2)
    */
    inline Double_t abs2() const { return re_ * re_ + im_ * im_; }

    //! Obtain the phase angle of the complex number
    /*!
        \return the phase angle of the complex number
    */
    inline Double_t arg() const { return TMath::ATan2( im_, re_ ); }

    //! Obtain the exponential of the complex number
    /*!
        \return the exponential of the complex number
    */
    inline LauComplex exp() const
    {
        Double_t mag( TMath::Exp( re_ ) );
        return LauComplex( mag * TMath::Cos( im_ ), mag * TMath::Sin( im_ ) );
    }

    //! Obtain the complex conjugate
    /*!
        \return the complex conjugate
    */
    inline LauComplex conj() const { return LauComplex( re_, -im_ ); }

    //! Transform this to its complex conjugate
    inline void makeConj() { im_ = -im_; }

    //! Obtain the complex number scaled by some factor
    /*!
        \param [in] scaleVal the value used to scale the complex number
        \return the complex number scaled by the scaleVal
    */
    inline LauComplex scale( Double_t scaleVal ) const
    {
        return LauComplex( scaleVal * re_, scaleVal * im_ );
    }

    //! Scale this by a factor
    /*!
        \param [in] scaleVal the value used to scale the complex number
    */
    inline void rescale( Double_t scaleVal )
    {
        re_ *= scaleVal;
        im_ *= scaleVal;
    }

    //! Set the real part
    /*!
        \param [in] realpart the value to be set as the real part
    */
    inline void setRealPart( Double_t realpart ) { re_ = realpart; }

    //! Set the imaginary part
    /*!
        \param [in] imagpart the value to be set as the imaginary part
    */
    inline void setImagPart( Double_t imagpart ) { im_ = imagpart; }

    //! Set both real and imaginary part
    /*!
        \param [in] realpart the value to be set as the real part
        \param [in] imagpart the value to be set as the imaginary part
    */
    inline void setRealImagPart( Double_t realpart, Double_t imagpart )
    {
        this->setRealPart( realpart );
        this->setImagPart( imagpart );
    }

    //! Set both real and imaginary part to zero
    inline void zero() { this->setRealImagPart( 0.0, 0.0 ); }

    //! Print the complex number
    void print() const;

  private:
    //! The real part
    Double_t re_;
    //! The imaginary part
    Double_t im_;

    ClassDef( LauComplex, 0 )    // a non-persistent bare-bones complex class
};

//! input operator formatting of a complex number
std::istream& operator>>( std::istream& os, LauComplex& z );
//! output operator formatting of a complex number
std::ostream& operator<<( std::ostream& os, const LauComplex& z );

#endif
