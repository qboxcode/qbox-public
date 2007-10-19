////////////////////////////////////////////////////////////////////////////////
//
// LDAFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LDAFunctional.h,v 1.5 2007-10-19 16:24:04 fgygi Exp $

#ifndef LDAFUNCTIONAL_H
#define LDAFUNCTIONAL_H

#include <vector>
#include <cassert>
#include "XCFunctional.h"

class LDAFunctional : public XCFunctional
{
  void xc_unpolarized(const double rh, double &ee, double &vv);
  void xc_polarized(const double rh, double &ee, double &vv);
  std::vector<double> _exc;
  std::vector<std::vector<double> > _vxc;

  LDAFunctional();

  public:

  LDAFunctional(const std::vector<std::vector<double> > &rhoe)
  {
    _nspin = rhoe.size();
    if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
    _np = rhoe[0].size();
    _exc.resize(_np);
    _vxc.resize(_nspin);
    for ( int i = 0; i < _nspin; i++ )
    {
      _vxc[i].resize(_np);
    }

    if ( _nspin == 1 )
    {
      rho = &rhoe[0][0];
      exc = &_exc[0];
      vxc1 = &_vxc[0][0];
    }
    else
    {
      rho_up = &rhoe[0][0];
      rho_dn = &rhoe[1][0];
      exc = &_exc[0];
      vxc1_up = &_vxc[0][0];
      vxc1_dn = &_vxc[1][0];
    }
  };

  bool isGGA() { return false; };
  std::string name() { return "LDA"; };
  void setxc(void);
};
#endif
