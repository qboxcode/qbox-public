////////////////////////////////////////////////////////////////////////////////
//
// Control.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Control.h,v 1.5 2003-12-02 20:24:27 fgygi Exp $

#ifndef CONTROL_H
#define CONTROL_H

#include <string>

struct Control
{
  // control variables
  string wf_dyn, atoms_dyn, cell_dyn; // dynamics string flags 
  int nite;
  string fermi;        // use Fermi distribution to fill states
  double ecutprec;
  double ecuts,sigmas,facs; // confinement energy parameters
  double prefmbar;      // externally applied pressure (Mbar)

  string wf_diag;
  string lock_abc;
  
  string tcp;
  double tcp_rcut;
  double tcp_sigma;
  
  double gms_mix; // mixing factor for generalized minimum spread functions
  
  string solvation; // continuum dielectric model for solvent

  string thermostat;
  string stress;
  string xc;
  string spin;
  int delta_spin;

  double dt;
  int iprint;
  int timeout; 

  double th_temp,th_time; // thermostat control
  double fermi_temp;  // temperature of Fermi distribution
  double emass;       // electron mass
  double vmass;       // cell mass
  
};
#endif
