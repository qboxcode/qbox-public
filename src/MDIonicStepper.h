////////////////////////////////////////////////////////////////////////////////
//
// MDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDIonicStepper.h,v 1.4 2004-03-11 21:52:31 fgygi Exp $

//
// IonicStepper is used in the following way
//
// input: r0,v0
//
// for ( k=0; k<n; k++ )
// {
//   compute forces f0(r0)
//   if ( k > 0 ) stepper->compute_v0(f0);
//   optional: modify velocities
//   stepper->update_v();
//   // new r0,v0,f0 known
//   stepper->compute_rp(f0); // using r0, v0 and f0
//   optional: restore constraints using rp and r0
//   stepper->update_r(); // r0 <- rp
//   // Note: r0 and v0 are not consistent at this point
// }
// compute forces f0(r0)
// stepper->compute_v0(f0);
// stepper->update_v();
//
// // r0,v0,f0 consistent at this point
//

#ifndef MDIONICSTEPPER_H
#define MDIONICSTEPPER_H

#include "IonicStepper.h"

class MDIonicStepper : public IonicStepper
{
  private:
  
  double th_temp_;
  double th_time_;
  double th_width_;
  double eta_;
  bool thermostat_;
  vector<vector< double> >  vhalf_;      // vhalf_[nsp_][3*na_]: v(t+dt/2)

  public:
  
  MDIonicStepper(Sample& s) : IonicStepper(s)
  {
    thermostat_ = ( s.ctrl.thermostat == "ON" );
    th_temp_ = s.ctrl.th_temp;
    th_time_ = s.ctrl.th_time;
    th_width_ = s.ctrl.th_width;
    eta_ = 0.0;
    vhalf_.resize(nsp_);
    for ( int is = 0; is < nsp_; is++ )
    {
      const int nais = atoms_.na(is);
      vhalf_[is].resize(3*nais);
    }
    atoms_.get_velocities(v0_);
  }

  void compute_rp(const vector<vector< double> >& f0);
  void compute_v0(const vector<vector< double> >& f0);
  void update_r(void);
  void update_v(void);
  double eta(void) const { return eta_; }
};

#endif
