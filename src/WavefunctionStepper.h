////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionStepper.h,v 1.2 2003-11-21 20:01:47 fgygi Exp $

#ifndef WAVEFUNCTIONSTEPPER_H
#define WAVEFUNCTIONSTEPPER_H
class Wavefunction;

class WavefunctionStepper
{
  private:
  
  public:

  virtual void update(Wavefunction& dwf) = 0;

  WavefunctionStepper() {}
  virtual ~WavefunctionStepper() {}
};
#endif
