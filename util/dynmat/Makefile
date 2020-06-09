#-------------------------------------------------------------------------------
#
# Copyright (c) 2008-2014 The Regents of the University of California
#
# This file is part of Qbox
#
# Qbox is distributed under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 2 of
# the License, or (at your option) any later version.
# See the file COPYING in the root directory of this distribution
# or <http://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------
#
EXEC=dynmat
LDFLAGS=-llapack
CXXFLAGS=-g
$(EXEC): dynmat.cpp
	$(CXX) -o $@ $^ $(LDFLAGS)
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f *.o
	rm -f $(EXEC)
#------------------------------------------------------------------------------
