////////////////////////////////////////////////////////////////////////////////
//
// Control.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Control.h,v 1.6 2004-02-04 19:55:17 fgygi Exp $

#ifndef CONTROL_H
#define CONTROL_H

#include <string>
#include <vector>

struct Control
{
  // control variables
  string wf_dyn, atoms_dyn; // dynamics string flags 
  int nite;
  double emass;       // electron mass
  
  string fermi;        // use Fermi distribution to fill states
  double fermi_temp;  // temperature of Fermi distribution
  double ecutprec;

  string wf_diag;
  
  string tcp;
  double tcp_rcut;
  double tcp_sigma;
  
  double gms_mix; // mixing factor for generalized minimum spread functions
  
  string thermostat;
  double th_temp,th_time; // thermostat control
  
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
