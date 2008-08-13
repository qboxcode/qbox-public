////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 2 of 
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// XMLGFPreprocessor.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: XMLGFPreprocessor.h,v 1.4 2008-08-13 06:39:43 fgygi Exp $

#include <string>
#include "Matrix.h"

////////////////////////////////////////////////////////////////////////////////
//
// XMLGFPreprocessor class
// Used to preprocess <grid_function> elements in an XML document
// Input: filename, DoubleMatrix& gfdata, string xmlcontent
// The preprocessor reads the file in parallel, processes all <grid_function>
// elements and stores the values of the grid_functions in the matrix gfdata
// which has dimensions (ngf,maxgridsize), where ngf is the total number of
// <grid_function> elements found in the file, and maxgridsize is the size of
// the largest grid_function.
// On return, the string xmlcontent contains (on task 0) the XML file
// with <grid_function> elements reduced to empty strings.
//
////////////////////////////////////////////////////////////////////////////////
class XMLGFPreprocessor
{
  public:

  void process(const char* const filename,
    DoubleMatrix& gfdata, std::string& xmlcontent);
};
