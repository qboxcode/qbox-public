////////////////////////////////////////////////////////////////////////////////
//
// StrX.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: StrX.h,v 1.2 2003-05-16 16:14:00 fgygi Exp $

#ifndef STRX_H
#define STRX_H

#include <string>
#include <iostream>
#include <cstdlib>
using namespace std;
#include <xercesc/util/XMLString.hpp>
using namespace xercesc;

class StrX
{
public :
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    StrX(const XMLCh* const toTranscode)
    {
        // Call the private transcoding method
        fLocalForm = XMLString::transcode(toTranscode);
    }

    ~StrX()
    {
        delete [] fLocalForm;
    }

    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const char* localForm() const
    {
        return fLocalForm;
    }

private :
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fLocalForm
    //      This is the local code page form of the string.
    // -----------------------------------------------------------------------
    char*   fLocalForm;
};

inline ostream& operator<<(ostream& target, const StrX& toDump)
{
    target << toDump.localForm();
    return target;
}
#endif
