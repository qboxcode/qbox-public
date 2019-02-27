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
///////////////////////////////////////////////////////////////////////////////
//
// readTSC.C
//
///////////////////////////////////////////////////////////////////////////////

unsigned long long readTSC(void)
{
  union { unsigned long long complete; unsigned int part[2]; } ticks;
#ifdef __amd64__
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
            : "=mr" (ticks.part[0]),
              "=mr" (ticks.part[1])
            : /* no inputs */
            : "eax", "edx");
#elif __powerpc__
  unsigned int tmp;
  __asm__ ("0:"
            "mftbu %[hi32]\n"
            "mftb %[lo32]\n"
            "mftbu %[tmp]\n"
            "cmpw %[tmp],%[hi32]\n"
            "bne 0b\n"
            : [hi32] "=r"(ticks.part[0]), [lo32] "=r"(ticks.part[1]),
            [tmp] "=r"(tmp));
#endif
  return ticks.complete;
}
