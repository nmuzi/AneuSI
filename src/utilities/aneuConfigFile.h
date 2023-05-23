/*=========================================================================
Module    : AneuSurfIsolation
File      : aneuConfigFile.h
Copyright : (C)opyright 2009++
            See COPYRIGHT statement in top level directory.
Authors   : Richard J. Wagner
Modified  : Daniel Millan
Purpose   : Class for reading named values from configuration files
Date      :
Version   :
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
/** 
\brief Class for reading named values from configuration files
Richard J. Wagner  v2.1  24 May 2004  wagnerr@umich.edu

Copyright (c) 2004 Richard J. Wagner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

Typical usage
-------------

Given a configuration file "settings.inp":
atoms  = 25
length = 8.0  # nanometers
name = Reece Surcher

Named values are read in various ways, with or without default values:
aneuConfigFile config( "settings.inp" );
int atoms = config.read<int>( "atoms" );
double length = config.read( "length", 10.0 );
string author, title;
config.readInto( author, "name" );
config.readInto( title, "title", string("Untitled") );

See file example.cpp for more examples.
*/
#ifndef __bioConfigfile_h
#define __bioConfigfile_h

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using std::string;

struct file_not_Found
{
    string filename;
};
typedef file_not_Found File_Not_Found;
    
struct key_not_Found
{   // thrown only by T read(key) variant of read()
    string key;
};
typedef key_not_Found Key_Not_Found;

    
class aneuConfigFile {
    // Data
    protected:
        string myDelimiter;  // separator between key and value
        string myComment;    // separator between value and comments
        string mySentry;     // optional string to signal end of file
        std::map<string,string> myContents;  // extracted keys and values
        
        typedef std::map<string,string>::iterator mapi;
        typedef std::map<string,string>::const_iterator mapci;

    // Methods
    public:
        aneuConfigFile( string filename,
                    string delimiter = "=",
                    string comment = "#",
                    string sentry = "EndConfigFile" );
        aneuConfigFile();
        
        // Search for key and read value or optional default value
        template<class T> T read( const string& key ) const;  // call as read<T>
        template<class T> T read( const string& key, const T& value ) const;
        template<class T> bool readInto( T& var, const string& key ) const;
        template<class T>
        bool readInto( T& var, const string& key, const T& value ) const;
        
        // Modify keys and values
        template<class T> void add( string key, const T& value );
        void remove( const string& key );
        
        // Check whether key exists in configuration
        bool keyExists( const string& key ) const;
        
        // Check or change configuration syntax
        string getDelimiter() const { return myDelimiter; }
        string getComment() const { return myComment; }
        string getSentry() const { return mySentry; }
        string setDelimiter( const string& s )
        { string old = myDelimiter;  myDelimiter = s;  return old; }  
        string setComment( const string& s )
        { string old = myComment;  myComment = s;  return old; }

        // Write or read configuration
        friend std::ostream& operator<<( std::ostream& os, const aneuConfigFile& cf );
        friend std::istream& operator>>( std::istream& is, aneuConfigFile& cf );

    protected:
        template<class T> static string T_as_string( const T& t );
        template<class T> static T string_as_T( const string& s );
        static void trim( string& s );


    // Exception types
    public:
//         File_Not_Found file_not_found( const string& filename_ )
//         {
//             File_Not_Found fnf;
//             fnf.filename = filename_;
//             return fnf;
//         }
//         File_Not_Found file_not_found(  )
//         {
//             File_Not_Found fnf;
//             fnf.filename =  string();
//             return fnf;
//         }
//         Key_Not_Found key_not_found( const string& key_ = string() )
//         {
//             Key_Not_Found knf;
//             knf.key = key_;
//             return knf;
//         }
        struct file_not_found {
            string filename;
            file_not_found( const string& filename_ = string() )
                : filename(filename_) {} };
        struct key_not_found {  // thrown only by T read(key) variant of read()
            string key;
            key_not_found( const string& key_ = string() )
                : key(key_) {} };
};


// static
template<class T>
string aneuConfigFile::T_as_string( const T& t )
{
    // Convert from a T to a string
    // Type T must support << operator
    std::ostringstream ost;
    ost << t;
    return ost.str();
}


/* static */
template<class T>
T aneuConfigFile::string_as_T( const string& s )
{
    // Convert from a string to a T
    // Type T must support >> operator
    T t;
    std::istringstream ist(s);
    ist >> t;
    return t;
}


/* static */
template<>
inline string aneuConfigFile::string_as_T<string>( const string& s )
{
    // Convert from a string to a string
    // In other words, do nothing
    return s;
}


/* static */
template<>
inline bool aneuConfigFile::string_as_T<bool>( const string& s )
{
    // Convert from a string to a bool
    // Interpret "false", "F", "no", "n", "0" as false
    // Interpret "true", "T", "yes", "y", "1", "-1", or anything else as true
    bool b = true;
    string sup = s;
    for( string::iterator p = sup.begin(); p != sup.end(); ++p )
        *p = (char) std::toupper(*p);  // The toupper function returns the uppercase version of the character
    if( sup==string("FALSE") || sup==string("F") ||
        sup==string("NO") || sup==string("N") ||
        sup==string("0") || sup==string("NONE") )
        b = false;
    return b;
}


template<class T>
T aneuConfigFile::read( const string& key ) const
{
    // Read the value corresponding to key
    mapci p = myContents.find(key);
    if( p == myContents.end() ) throw key_not_found(key);
    return string_as_T<T>( p->second );
}


template<class T>
T aneuConfigFile::read( const string& key, const T& value ) const
{
    // Return the value corresponding to key or given default value
    // if key is not found
    mapci p = myContents.find(key);
    if( p == myContents.end() ) return value;
    return string_as_T<T>( p->second );
}


template<class T>
bool aneuConfigFile::readInto( T& var, const string& key ) const
{
    // Get the value corresponding to key and store in var
    // Return true if key is found
    // Otherwise leave var untouched
    mapci p = myContents.find(key);
    bool found = ( p != myContents.end() );
    if( found ) var = string_as_T<T>( p->second );
    return found;
}


template<class T>
bool aneuConfigFile::readInto( T& var, const string& key, const T& value ) const
{
    // Get the value corresponding to key and store in var
    // Return true if key is found
    // Otherwise set var to given default
    mapci p = myContents.find(key);
    bool found = ( p != myContents.end() );
    if( found )
        var = string_as_T<T>( p->second );
    else
        var = value;
    return found;
}


template<class T>
void aneuConfigFile::add( string key, const T& value )
{
    // Add a key with given value
    string v = T_as_string( value );
    trim(key);
    trim(v);
    myContents[key] = v;
    return;
}

#endif  // CONFIGFILE_H

// Release notes:
// v1.0  21 May 1999
//   + First release
//   + Template read() access only through non-member readConfigFile()
//   + ConfigurationFileBool is only built-in helper class
// 
// v2.0  3 May 2002
//   + Shortened name from ConfigurationFile to aneuConfigFile
//   + Implemented template member functions
//   + Changed default comment separator from % to #
//   + Enabled reading of multiple-line values
// 
// v2.1  24 May 2004
//   + Made template specializations inline to avoid compiler-dependent linkage
//   + Allowed comments within multiple-line values
//   + Enabled blank line termination for multiple-line values
//   + Added optional sentry to detect end of configuration file
//   + Rewrote messy trimWhitespace() function as elegant trim()

// **************************************************************************************

// README for aneuConfigFile distribution
// Richard J. Wagner  v2.1  24 May 2004
// 
// Instructions
// ------------
// 
// The only necessary files for using this configuration file reader are
// "aneuConfigFile.h" and "aneuConfigFile.cpp".  The class name is aneuConfigFile.
// 
// Usage examples are in "example.cpp".  Linux or Unix users can type "make" to
// compile and then type "make run" to run the example program.
// 
// The test program in "tester.cpp" will check that the class properly reads
// a variety of simple and complex configuration file entries.  To run the test
// program type "make test".
// 
// When you are done with the examples and the test program, type "make clean"
// to get rid of temporary files.
// 
// For Windows or Mac users with a compiler such as Metrowerks CodeWarrior or
// Microsoft Visual C++, simply add "example.cpp" and "aneuConfigFile.cpp" to an
// empty C++ console application.  Compile and run to see the configuration
// file reader in action.  Do likewise with "tester.cpp" to check that the
// code works properly with your compiler.
// 
// If you encounter any problems, please e-mail a copy of the output and a
// description of the test system to me at "wagnerr@umich.edu".  Any other
// feedback is welcome too.
// 
// 
// Installation
// ------------
// 
// Just copy the files "aneuConfigFile.h" and "aneuConfigFile.cpp" to your working
// directory or some other place where your compiler can find them.  Add
// "aneuConfigFile.cpp" to your project and put the following line at the top of
// your program to access the aneuConfigFile class:
// 
// #include "aneuConfigFile.h"
// 
// 
// Contents
// --------
// 
// README            - this file
// aneuConfigFile.h      - declaration of aneuConfigFile class
// aneuConfigFile.cpp    - definitions of aneuConfigFile class
// example.cpp       - examples of using aneuConfigFile
// tester.cpp        - tests aneuConfigFile class
// example.inp       - configuration file for example program
// test.inp          - configuration file for tester program
// Triplet.h         - sample user-defined data type
// Makefile          - instructions used by "make" command
// aneuConfigFile.html   - Web page about aneuConfigFile
// AntBlueMaize.jpg  - background for aneuConfigFile.html
// ArrowHome.gif     - home icon for aneuConfigFile.html
// main.css          - style sheet for aneuConfigFile.html
