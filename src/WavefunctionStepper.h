////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionStepper.h,v 1.3 2003-11-27 01:20:59 fgygi Exp $

#ifndef WAVEFUNCTIONSTEPPER_H
#define WAVEFUNCTIONSTEPPER_H
class Wavefunction;

class WavefunctionStepper
{
  private:
  
  public:

  virtual void update(Wavefunction& dwf) = 0;
  virtual void preprocess(void) {}
  virtual void postprocess(void) {}

  WavefunctionStepper() {}
  virtual ~WavefunctionStepper() {}
};
#endif
