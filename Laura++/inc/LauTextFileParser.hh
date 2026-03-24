
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

/*! \file LauTextFileParser.hh
    \brief File containing declaration of LauTextFileParser class.
*/

/*! \class LauTextFileParser
    \brief Class for parsing text files.

    Class for parsing text files. Each line is split on spaces.
*/

#ifndef LAU_TEXT_FILE_PARSER_HH
#define LAU_TEXT_FILE_PARSER_HH

#include "TString.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class LauTextFileParser {

  public:
    //! Constructor
    /*!
        \param [in] s filename
    */
    LauTextFileParser( const char* s );

    //! Constructor
    /*!
        \param [in] s filename
    */
    LauTextFileParser( const std::string& s );

    //! Constructor
    /*!
        \param [in] s filename
    */
    LauTextFileParser( const TString& s );

    //! Destructor
    virtual ~LauTextFileParser();

    //! Parse the file
    void processFile();

    //! Retrieve the specified line
    /*!
        \param [in] lineNo the line number
        \return the specified line
    */
    std::vector<std::string> getLine( UInt_t lineNo );

    //! Retrieve the next line
    /*!
        \return the next line
    */
    std::vector<std::string> getNextLine();

    //! Set the line number
    /*!
        \param [in] lineNo the new line number
    */
    void resetLineNumber( UInt_t lineNo ) { lineNumber_ = lineNo; }

    //! Get the total number of lines that are not comments
    /*!
        \return the total number of non-comment lines
    */
    UInt_t getTotalNumLines() const { return totalLines_; }

  protected:

  private:
    //! Copy constructor (not implemented)
    LauTextFileParser( const LauTextFileParser& rhs );

    //! Copy assignment operator (not implemented)
    LauTextFileParser& operator=( const LauTextFileParser& rhs );

    //! The type used to contain the parsed text
    typedef std::map<int, std::vector<std::string>> LineMap;

    //! The name of the file
    std::string theFile_;

    //! The parsed text
    LineMap lineMap_;

    //! The total number of lines in the file
    UInt_t totalLines_;

    //! The current line number
    UInt_t lineNumber_;
};

#endif
