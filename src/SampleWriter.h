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
// SampleWriter.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SAMPLEWRITER_H
#define SAMPLEWRITER_H

#include "Context.h"
#include <string>
class Sample;

class SampleWriter
{
  public:

  void writeSample(const Sample& s, const std::string filename,
                   std::string description,
                   bool base64, bool atomsonly, bool serial, bool save_wfv);
};

class SampleWriterException
{
  public:
  std::string msg;
  SampleWriterException(std::string s) : msg(s) {}
};

#endif
