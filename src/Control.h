////////////////////////////////////////////////////////////////////////////////
//
// Control.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Control.h,v 1.8 2004-09-14 22:24:11 fgygi Exp $

#ifndef CONTROL_H
#define CONTROL_H

#include <string>
#include <vector>

struct Control
{
  // control variables
  string debug, timing;
  string wf_dyn, atoms_dyn; // dynamics string flags 
  int nite;
  double emass;       // electron mass
  
  double fermi_temp;  // temperature of Fermi distribution
  double ecutprec;

  string wf_diag;
  
  string tcp;
  double tcp_rcut;
  double tcp_sigma;
  
  double gms_mix; // mixing factor for generalized minimum spread functions
  
  string thermostat;
  double th_temp,th_time, th_width; // thermostat control
  
  string stress;
  string cell_dyn;
  string cell_lock;
  double cell_mass;
  double ecuts,sigmas,facs; // confinement energy parameters
  double ext_stress[6]; // external stress tensor: xx,yy,zz,xy,yz,xz
  
  string xc;
  string spin;
  int delta_spin;

  double dt;
  int iprint;
  int timeout; 
};
#endif
