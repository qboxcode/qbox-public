////////////////////////////////////////////////////////////////////////////////
//
// MDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDWavefunctionStepper.h,v 1.3 2004-03-11 21:52:31 fgygi Exp $

#ifndef MDWAVEFUNCTIONSTEPPER_H
#define MDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class MDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double dt_, dt2bye_;

  double ekin_ep_, ekin_em_;
  double ekin_eh(void);

  public:

  void update(Wavefunction& dwf);
  void compute_wfm(Wavefunction& dwf);
  void compute_wfv(Wavefunction& dwf);
  double ekin(void) const { return 0.5*(ekin_ep_ + ekin_em_); }

  MDWavefunctionStepper(Sample& s, TimerMap& tmap);
  ~MDWavefunctionStepper() {};
};
#endif
