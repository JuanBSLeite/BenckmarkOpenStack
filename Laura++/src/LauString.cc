
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

/*! \file LauString.cc
    \brief File containing implementation of LauString class.
*/

#include "LauString.hh"

#include <iostream>
#include <sstream>

using std::cout;
using std::endl;

LauString::LauString() :
    theString_( "" )
{
}

LauString::LauString( const char* s ) :
    theString_( s )
{
}

LauString::LauString( const std::string& s ) :
    theString_( s )
{
}

LauString::LauString( const TString& s ) :
    theString_( s.Data() )
{
}

LauString::LauString( const LauString& s ) :
    theString_( s.theString_ )
{
}

LauString& LauString::operator=( const LauString& s )
{
    if ( this != &s ) {
        theString_ = s.theString_;
    }
    return *this;
}

LauString::~LauString()
{
}

void LauString::addInteger( int integer )
{
    std::ostringstream o;
    o << integer;
    theString_ += o.str();
}

void LauString::addDouble( double value )
{
    std::ostringstream o;
    o << value;
    theString_ += o.str();
}

LauString LauString::operator+=( int integer )
{
    this->addInteger( integer );
    return ( *this );
}

LauString LauString::operator+=( double value )
{
    this->addDouble( value );
    return ( *this );
}

LauString LauString::operator+=( const char* text )
{
    this->addText( text );
    return ( *this );
}

LauString LauString::operator+=( const std::string& text )
{
    this->addText( text );
    return ( *this );
}

LauString LauString::operator+=( const TString& text )
{
    this->addText( text.Data() );
    return ( *this );
}

std::vector<std::string> LauString::split( const LauString& splitter ) const
{
    std::string splitString         = splitter.getString();
    std::vector<std::string> result = this->split( splitString );
    return result;
}

std::vector<std::string> LauString::split( const std::string& splitter ) const
{
    // Code from STLplus
    std::vector<std::string> result;

    if ( ! theString_.empty() && ! splitter.empty() ) {

        for ( std::string::size_type offset = 0;; ) {

            std::string::size_type found = theString_.find( splitter, offset );

            if ( found != std::string::npos ) {
                std::string tmpString = theString_.substr( offset, found - offset );
                if ( tmpString.size() > 0 ) {
                    result.push_back( tmpString );
                }
                offset = found + splitter.size();
            } else {
                std::string tmpString = theString_.substr( offset, theString_.size() - offset );
                if ( tmpString.size() > 0 ) {
                    result.push_back( tmpString );
                }
                break;
            }
        }
    }

    return result;
}
