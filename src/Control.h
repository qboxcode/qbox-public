////////////////////////////////////////////////////////////////////////////////
//
// Control.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Control.h,v 1.11 2007-03-17 01:14:00 fgygi Exp $

#ifndef CONTROL_H
#define CONTROL_H

#include <string>
#include <vector>

struct Control
{
  // control variables
  std::string debug, timing;
  std::string wf_dyn, atoms_dyn; // dynamics string flags 
  int nite;
  double emass;       // electron mass
  
  double fermi_temp;  // temperature of Fermi distribution
  double ecutprec;

  std::string wf_diag;
  
  std::string tcp;
  double tcp_rcut;
  double tcp_sigma;
  
  double gms_mix; // mixing factor for generalized minimum spread functions
  
  std::string thermostat;
  double th_temp,th_time, th_width; // thermostat control
  
  std::string stress;
  std::string cell_dyn;
  std::string cell_lock;
  double cell_mass;
  double ecuts;         // confinement potential energy cutoff
  double ext_stress[6]; // external stress tensor: xx,yy,zz,xy,yz,xz
  
  std::string xc;
  std::string spin;
  int delta_spin;

  double dt;
  int iprint;
  int timeout;
  
  double charge_mix_coeff;
  double charge_mix_rcut;
};
#endif
