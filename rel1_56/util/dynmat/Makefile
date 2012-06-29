#-------------------------------------------------------------------------------
#
# Copyright (c) 2008 The Regents of the University of California
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
# $Id: Makefile,v 1.65 2009/09/08 05:40:50 fgygi Exp $
#------------------------------------------------------------------------------
#
include $(TARGET).mk
#------------------------------------------------------------------------------
SRCDIR = $(HOME)/qb/export/rel1_52_1/qb/src/
EXEC = dynmat
# use source files in SRCDIR but not objects, which may use MPI
vpath %.h $(SRCDIR)
vpath %.C $(SRCDIR)
OBJECTS = dynmat.o Matrix.o Context.o
CXXFLAGS += -DTARGET='"$(TARGET)"' -I$(SRCDIR)
 $(EXEC):    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
 depend :
	makedepend -Y -D$(PLT) *.[cCh]
#makedepend -D$(PLT) -I/usr/include/CC *.[cC]
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f *.o
	rm -f $(EXEC) $(TESTEXECS)
#------------------------------------------------------------------------------
 ctags :
	etags -o tags *.[Ch]
#------------------------------------------------------------------------------
# DO NOT DELETE

Context.o: Context.h
dynmat.o: Context.h Matrix.h
Matrix.o: Context.h Matrix.h
