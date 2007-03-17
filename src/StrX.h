////////////////////////////////////////////////////////////////////////////////
//
// StrX.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: StrX.h,v 1.3 2007-03-17 01:14:00 fgygi Exp $

#ifndef STRX_H
#define STRX_H

#include <string>
#include <iostream>
#include <cstdlib>
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

inline std::ostream& operator<<(std::ostream& target, const StrX& toDump)
{
    target << toDump.localForm();
    return target;
}
#endif
