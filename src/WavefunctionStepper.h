////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionStepper.h,v 1.1.1.1 2002-09-27 00:08:39 fgygi Exp $

#ifndef WAVEFUNCTIONSTEPPER_H
#define WAVEFUNCTIONSTEPPER_H

class WavefunctionStepper
{
  private:

  public:

  virtual void update(Wavefunction& wf, WavefunctionDerivative& dwf) = 0;

  WavefunctionStepper();
  ~WavefunctionStepper();
};
#endif
