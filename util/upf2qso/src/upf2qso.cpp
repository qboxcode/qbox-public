////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009-2017 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// upf2qso.cpp: transform a UPF pseudopotential to QSO format
// the QSO format is defined at http://www.quantum-simulation.org
// This program accepts both the original UPF format and the "2.0.1" format
//
////////////////////////////////////////////////////////////////////////////////
//
// use: upf2qso rcut(a.u.) < pot.UPF > pot.xml
// The potentials in pot.xml are given on a radial linear grid up to rcut,
// with grid spacing of 0.01 (a.u.)
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
#include "spline.h"
#include "isodate.h"
#include "PeriodicTable.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
string find_start_element(string name)
{
  // return the contents of the tag at start of element "name"
  string buf, token;
  string search_str = "<" + name;
  do
  {
    cin >> token;
  }
  while ( !cin.eof() && token.find(search_str) == string::npos );
  if ( cin.eof() )
  {
    cerr << " EOF reached before start element " << name << endl;
    throw invalid_argument(name);
  }

  buf = token;
  if ( buf[buf.length()-1] == '>' )
    return buf;

  // read until ">" is found
  bool found = false;
  char ch;
  do
  {
    cin.get(ch);
    found = ch == '>';
    buf += ch;
  }
  while ( !cin.eof() && !found );
  if ( cin.eof() )
  {
    cerr << " EOF reached before > " << name << endl;
    throw invalid_argument(name);
  }
  return buf;
}

////////////////////////////////////////////////////////////////////////////////
void find_end_element(string name)
{
  string buf, token;
  string search_str = "</" + name + ">";
  do
  {
    cin >> token;
    if ( token.find(search_str) != string::npos ) return;
  }
  while ( !cin.eof() );
  cerr << " EOF reached before end element " << name << endl;
  throw invalid_argument(name);
}

////////////////////////////////////////////////////////////////////////////////
void seek_str(string tag)
{
  // Read tokens from stdin until tag is found.
  // Throw an exception if tag not found before eof()
  bool done = false;
  string token;
  int count = 0;

  do
  {
    cin >> token;
    if ( token.find(tag) != string::npos ) return;
  }
  while ( !cin.eof() );

  cerr << " EOF reached before " << tag << endl;
  throw invalid_argument(tag);
}

////////////////////////////////////////////////////////////////////////////////
string get_attr(string buf, string attr)
{
  //cerr << " get_attr: searching for: " << attr << endl;
  //cerr << " in buffer: " << endl;
  //cerr << buf << endl;
  bool done = false;
  string s, search_string = attr + "=";

  // find attribute name in buf
  string::size_type p = buf.find(search_string);
  if ( p != string::npos )
  {
    // process attribute
    string::size_type b = buf.find_first_of("\"",p);
    string::size_type e = buf.find_first_of("\"",b+1);
    if ( b == string::npos || e == string::npos )
    {
      cerr << " get_attr: attribute not found: " << attr << endl;
      throw invalid_argument(attr);
    }
    s = buf.substr(b+1,e-b-1);
    // remove white space
    s.erase(remove_if(s.begin(),s.end(),::isspace),s.end());
    return s;
  }
  else
  {
    cerr << " get_attr: attribute not found: " << attr << endl;
    throw invalid_argument(attr);
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////////
void skipln(void)
{
  char ch;
  bool found = false;
  while ( !cin.eof() && !found )
  {
    cin.get(ch);
    found = ch == '\n';
  }
}

////////////////////////////////////////////////////////////////////////////////
const string release="1.10";

int main(int argc, char** argv)
{
  cerr << " upf2qso " << release << endl;
  if ( argc != 2 )
  {
    cerr << " usage: upf2qso rcut < file.upf > file.xml" << endl;
    return 1;
  }
  assert(argc==2);
  const double rcut = atof(argv[1]);

  PeriodicTable pt;
  string buf,s;
  istringstream is;

  // determine UPF version
  int upf_version = 0;

  // The first line of the UPF potential file contains either of the following:
  // <PP_INFO>  (for UPF version 1)
  // <UPF version="2.0.1"> (for UPF version 2)
  string::size_type p;
  getline(cin,buf);
  p = buf.find("<PP_INFO>");
  if ( p != string::npos )
    upf_version = 1;
  else
  {
    p = buf.find("<UPF version=\"2.0.1\">");
    if ( p != string::npos )
      upf_version = 2;
  }
  if ( upf_version == 0 )
  {
    cerr << " Format of UPF file not recognized " << endl;
    cerr << " First line of file: " << buf << endl;
    return 1;
  }

  cerr << " UPF version: " << upf_version << endl;

  if ( upf_version == 1 )
  {
    // process UPF version 1 potential
    string upf_pp_info;
    bool done = false;
    while (!done)
    {
      getline(cin,buf);
      is.clear();
      is.str(buf);
      is >> s;
      done = ( s == "</PP_INFO>" );
      if ( !done )
      {
        upf_pp_info += buf + '\n';
      }
    }

    // remove all '<' and '>' characters from the PP_INFO field
    // for XML compatibility
    p = upf_pp_info.find_first_of("<>");
    while ( p != string::npos )
    {
      upf_pp_info[p] = ' ';
      p = upf_pp_info.find_first_of("<>");
    }

    seek_str("<PP_HEADER>");
    skipln();

    // version number (ignore)
    getline(cin,buf);

    // element symbol
    string upf_symbol;
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_symbol;

    // get atomic number and mass
    const int atomic_number = pt.z(upf_symbol);
    const double mass = pt.mass(upf_symbol);

    // NC flag
    string upf_ncflag;
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_ncflag;
    if ( upf_ncflag != "NC" )
    {
      cerr << " not a Norm-conserving potential" << endl;
      cerr << " NC flag: " << upf_ncflag << endl;
      return 1;
    }

    // NLCC flag
    string upf_nlcc_flag;
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_nlcc_flag;
    if ( upf_nlcc_flag == "T" )
    {
      cerr << " Potential includes a non-linear core correction" << endl;
    }

    // XC functional (add in description)
    string upf_xcf[4];
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_xcf[0] >> upf_xcf[1] >> upf_xcf[2] >> upf_xcf[3];

    // add XC functional information to description
    upf_pp_info += upf_xcf[0] + ' ' + upf_xcf[1] + ' ' +
                   upf_xcf[2] + ' ' + upf_xcf[3] + '\n';

    // Z valence
    double upf_zval;
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_zval;

    // Total energy (ignore)
    getline(cin,buf);

    // suggested cutoff (ignore)
    getline(cin,buf);

    // max angular momentum
    int upf_lmax;
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_lmax;

    // number of points in mesh
    int upf_mesh_size;
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_mesh_size;

    // number of wavefunctions, number of projectors
    int upf_nwf, upf_nproj;
    getline(cin,buf);
    is.clear();
    is.str(buf);
    is >> upf_nwf >> upf_nproj;

    // Wavefunctions
    vector<string> upf_shell(upf_nwf);
    vector<int> upf_l(upf_nwf);
    vector<double> upf_occ(upf_nwf);
    // skip header
    getline(cin,buf);
    for ( int ip = 0; ip < upf_nwf; ip++ )
    {
      getline(cin,buf);
      is.clear();
      is.str(buf);
      is >> upf_shell[ip] >> upf_l[ip] >> upf_occ[ip];
    }
    seek_str("</PP_HEADER>");

    // read mesh
    seek_str("<PP_MESH>");
    seek_str("<PP_R>");
    skipln();
    vector<double> upf_r(upf_mesh_size);
    for ( int i = 0; i < upf_mesh_size; i++ )
     cin >> upf_r[i];
    seek_str("</PP_R>");
    seek_str("<PP_RAB>");
    skipln();
    vector<double> upf_rab(upf_mesh_size);
    for ( int i = 0; i < upf_mesh_size; i++ )
     cin >> upf_rab[i];
    seek_str("</PP_RAB>");
    seek_str("</PP_MESH>");

    vector<double> upf_nlcc;
    if ( upf_nlcc_flag == "T" )
    {
      upf_nlcc.resize(upf_mesh_size);
      seek_str("<PP_NLCC>");
      skipln();
      for ( int i = 0; i < upf_mesh_size; i++ )
        cin >> upf_nlcc[i];
      seek_str("</PP_NLCC>");
    }

    seek_str("<PP_LOCAL>");
    skipln();
    vector<double> upf_vloc(upf_mesh_size);
    for ( int i = 0; i < upf_mesh_size; i++ )
      cin >> upf_vloc[i];
    seek_str("</PP_LOCAL>");

    seek_str("<PP_NONLOCAL>");
    skipln();
    vector<vector<double> > upf_vnl;
    upf_vnl.resize(upf_nproj);
    vector<int> upf_proj_l(upf_nproj);
    for ( int j = 0; j < upf_nproj; j++ )
    {
      seek_str("<PP_BETA>");
      skipln();
      int ip, l, np;
      cin >> ip >> l;
      skipln();
      assert(ip-1 < upf_nproj);
      assert(l <= upf_lmax);
      upf_proj_l[ip-1] = l;
      cin >> np;
      upf_vnl[j].resize(upf_mesh_size);
      for ( int i = 0; i < np; i++ )
        cin >> upf_vnl[j][i];
      seek_str("</PP_BETA>");
      skipln();
    }
    seek_str("<PP_DIJ>");
    skipln();
    int upf_ndij;
    cin >> upf_ndij;
    skipln();
    if ( upf_ndij != upf_nproj )
    {
      cerr << " Number of non-zero Dij differs from number of projectors"
           << endl;
      return 1;
    }

    vector<double> upf_d(upf_ndij);
    for ( int i = 0; i < upf_ndij; i++ )
    {
      int m,n;
      cin >> m >> n >> upf_d[i];
      if ( m != n )
      {
        cerr << " Non-local Dij has off-diagonal elements" << endl;
        cerr << " m=" << m << " n=" << n << endl;
        return 1;
      }
    }
    seek_str("</PP_DIJ>");

    seek_str("</PP_NONLOCAL>");

    // make table iproj[l] mapping l to iproj
    // vnl(l) is in vnl[iproj[l]] if iproj[l] > -1
    // vlocal if iproj[llocal] = -1
    vector<int> iproj(upf_lmax+2);
    for ( int l = 0; l <= upf_lmax+1; l++ )
      iproj[l] = -1;
    for ( int j = 0; j < upf_nproj; j++ )
      iproj[upf_proj_l[j]] = j;

    // determine angular momentum of local potential in UPF file
    int upf_llocal;
    // reverse loop to get correct upf_llocal when upf_nproj < upf_lmax
    for ( int l = upf_lmax+1; l >= 0; l-- )
      if ( iproj[l] == -1 )
        upf_llocal = l;
    // increase lmax if there are more projectors than wavefunctions
    int qso_lmax = upf_lmax;
    if (upf_lmax < upf_llocal)
    {
      qso_lmax = upf_lmax+1;
    }

    seek_str("<PP_PSWFC>");
    skipln();
    vector<vector<double> > upf_wf;
    vector<int> upf_wf_l(upf_nwf);
    vector<double> upf_wf_occ(upf_nwf);
    upf_wf.resize(upf_nwf);
    for ( int j = 0; j < upf_nwf; j++ )
    {
      upf_wf[j].resize(upf_mesh_size);
      string label;
      cin >> label >> upf_wf_l[j] >> upf_wf_occ[j];
      skipln();
      for ( int i = 0; i < upf_mesh_size; i++ )
        cin >> upf_wf[j][i];
    }
    seek_str("</PP_PSWFC>");

    // output original data in file upf.dat
    ofstream upf("upf.dat");
    upf << "# vloc" << endl;
    for ( int i = 0; i < upf_vloc.size(); i++ )
      upf << upf_r[i] << " " << upf_vloc[i] << endl;
    upf << endl << endl;
    for ( int j = 0; j < upf_nproj; j++ )
    {
      upf << "# proj j=" << j << endl;
      for ( int i = 0; i < upf_vnl[j].size(); i++ )
        upf << upf_r[i] << " " << upf_vnl[j][i] << endl;
      upf << endl << endl;
    }
    for ( int j = 0; j < upf_nwf; j++ )
    {
      upf << "# wf j=" << j << endl;
      for ( int i = 0; i < upf_wf[j].size(); i++ )
        upf << upf_r[i] << " " << upf_wf[j][i] << endl;
      upf << endl << endl;
    }
    upf.close();


    // print summary
    cerr << "PP_INFO:" << endl << upf_pp_info << endl;
    cerr << "Element: " << upf_symbol << endl;
    cerr << "NC: " << upf_ncflag << endl;
    cerr << "NLCC: " << upf_nlcc_flag << endl;
    cerr << "XC: " << upf_xcf[0] << " " << upf_xcf[1] << " "
         << upf_xcf[2] << " " << upf_xcf[3] << endl;
    cerr << "Zv: " << upf_zval << endl;
    cerr << "lmax: " << qso_lmax << endl;
    cerr << "llocal: " << upf_llocal << endl;
    cerr << "nwf: " << upf_nwf << endl;
    cerr << "mesh_size: " << upf_mesh_size << endl;

    // compute delta_vnl[l][i] on the upf log mesh

    // divide the projector function by the wavefunction, except if
    // the wavefunction amplitude is smaller than tol, outside of rcut_divide.
    const double tol = 1.e-5;
    const double rcut_divide = 1.0;
    vector<vector<double> > delta_vnl;
    delta_vnl.resize(upf_nproj);
    for ( int j = 0; j < upf_nproj; j++ )
    {
      delta_vnl[j].resize(upf_wf[j].size());
      for ( int i = 0; i < delta_vnl[j].size(); i++ )
      {
        double num = upf_vnl[j][i];
        double den = upf_wf[upf_proj_l[j]][i];

        delta_vnl[j][i] = 0.0;
        if ( upf_r[i] < rcut_divide )
        {
          // near the origin
          if ( i == 0 && fabs(den) < tol )
          {
            // i = 0 for linear mesh, r = 0.0: use nearest value
            delta_vnl[j][i] = upf_vnl[j][1] / upf_wf[upf_proj_l[j]][1];
          }
          else
          {
            // other points near the origin
            delta_vnl[j][i] = num / den;
          }
        }
        else
        {
          // wavefunction gets small at large r.
          // Assume that delta_vnl is zero when that happens
          if ( fabs(den) > tol )
            delta_vnl[j][i] = num / den;
        }
      }
    }

    vector<vector<double> > vps;
    vps.resize(upf_nproj+1);
    for ( int j = 0; j < upf_nproj; j++ )
    {
      vps[j].resize(upf_mesh_size);
      for ( int i = 0; i < delta_vnl[j].size(); i++ )
        vps[j][i] = upf_vloc[i] + delta_vnl[j][i];
    }

    // interpolate functions on linear mesh
    const double mesh_spacing = 0.01;
    int nplin = (int) (rcut / mesh_spacing);
    vector<double> f(upf_mesh_size), fspl(upf_mesh_size);

    vector<double> nlcc_lin(nplin);
    // interpolate NLCC
    if ( upf_nlcc_flag == "T" )
    {
      assert(upf_mesh_size==upf_nlcc.size());
      for ( int i = 0; i < upf_nlcc.size(); i++ )
        f[i] = upf_nlcc[i];
      int n = upf_nlcc.size();
      int bcnat_left = 0;
      double yp_left = 0.0;
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);

      for ( int i = 0; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        if ( r >= upf_r[0] )
          splint(n,&upf_r[0],&f[0],&fspl[0],r,&nlcc_lin[i]);
        else
          // use value closest to the origin for r=0
          nlcc_lin[i] = upf_nlcc[0];
        if ( fabs(nlcc_lin[i]) < 1.e-12 )
          nlcc_lin[i] = 0.0;
      }
    }

    // interpolate vloc
    // factor 0.5: convert from Ry in UPF to Hartree in QSO
    for ( int i = 0; i < upf_vloc.size(); i++ )
      f[i] = 0.5 * upf_vloc[i];

    int n = upf_vloc.size();
    int bcnat_left = 0;
    double yp_left = 0.0;
    int bcnat_right = 1;
    double yp_right = 0.0;
    spline(n,&upf_r[0],&f[0],yp_left,yp_right,
           bcnat_left,bcnat_right,&fspl[0]);

    vector<double> vloc_lin(nplin);
    for ( int i = 0; i < nplin; i++ )
    {
      double r = i * mesh_spacing;
      if ( r >= upf_r[0] )
        splint(n,&upf_r[0],&f[0],&fspl[0],r,&vloc_lin[i]);
      else
        // use value closest to the origin for r=0
        vloc_lin[i] = 0.5 * upf_vloc[0];
    }

    // interpolate vps[j], j=0, nproj-1
    vector<vector<double> > vps_lin;
    vps_lin.resize(vps.size());
    for ( int j = 0; j < vps.size(); j++ )
    {
      vps_lin[j].resize(nplin);
    }

    for ( int j = 0; j < upf_nproj; j++ )
    {
      // factor 0.5: convert from Ry in UPF to Hartree in QSO
      for ( int i = 0; i < upf_vloc.size(); i++ )
        f[i] = 0.5 * vps[j][i];

      int n = upf_vloc.size();
      int bcnat_left = 0;
      double yp_left = 0.0;
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);

      for ( int i = 0; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        if ( r >= upf_r[0] )
          splint(n,&upf_r[0],&f[0],&fspl[0],r,&vps_lin[j][i]);
        else
          vps_lin[j][i] = 0.5 * vps[j][0];
      }
    }

    // write potentials in gnuplot format on file vlin.dat
    ofstream vlin("vlin.dat");
    for ( int l = 0; l <= qso_lmax; l++ )
    {
      vlin << "# v, l=" << l << endl;
      if ( iproj[l] == -1 )
      {
        // l == llocal
        for ( int i = 0; i < nplin; i++ )
          vlin << i*mesh_spacing << " " << vloc_lin[i] << endl;
        vlin << endl << endl;
      }
      else
      {
        for ( int i = 0; i < nplin; i++ )
          vlin << i*mesh_spacing << " " << vps_lin[iproj[l]][i] << endl;
        vlin << endl << endl;
      }
    }

    // interpolate wavefunctions on the linear mesh

    vector<vector<double> > wf_lin;
    wf_lin.resize(upf_nwf);
    for ( int j = 0; j < upf_nwf; j++ )
    {
      wf_lin[j].resize(nplin);
      assert(upf_wf[j].size()<=f.size());
      for ( int i = 0; i < upf_wf[j].size(); i++ )
      {
        if ( upf_r[i] > 0.0 )
          f[i] = upf_wf[j][i] / upf_r[i];
        else
        {
          // value at origin, depending on angular momentum
          if ( upf_wf_l[j] == 0 )
          {
            // l=0: take value closest to r=0
            f[i] = upf_wf[j][1]/upf_r[1];
          }
          else
          {
            // l>0:
            f[i] = 0.0;
          }
        }
      }

      int n = upf_wf[j].size();
      // choose boundary condition at origin depending on angular momentum
      int bcnat_left = 0;
      double yp_left = 0.0;
      if ( upf_wf_l[j] == 1 )
      {
        bcnat_left = 1; // use natural bc
        double yp_left = 0.0; // not used
      }
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);

      for ( int i = 0; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        if ( r >= upf_r[0] )
          splint(n,&upf_r[0],&f[0],&fspl[0],r,&wf_lin[j][i]);
        else
        {
          // r < upf_r[0]
          assert(upf_r[0]>0.0);
          // compute value near origin, depending on angular momentum
          if ( upf_wf_l[j] == 0 )
          {
            // l=0: take value closest to r=0
            wf_lin[j][i] = upf_wf[j][0]/upf_r[0];
          }
          else
          {
            // l>0:
            wf_lin[j][i] = upf_wf[j][0] * r / ( upf_r[0] * upf_r[0] );
          }
        }
      }

      vlin << "# phi, l=" << upf_l[j] << endl;
      for ( int i = 0; i < nplin; i++ )
        vlin << i*mesh_spacing << " " << wf_lin[j][i] << endl;
      vlin << endl << endl;
    }

    cerr << " interpolation done" << endl;

#if 1
    // output potential on log mesh
    ofstream vout("v.dat");
    for ( int l = 0; l <= qso_lmax; l++ )
    {
      vout << "# v, l=" << l << endl;
      if ( iproj[l] == -1 )
      {
        // l == llocal
        for ( int i = 0; i < upf_vloc.size(); i++ )
          vout << upf_r[i] << " " << 0.5*upf_vloc[i] << endl;
        vout << endl << endl;
      }
      else
      {
        for ( int i = 0; i < vps[iproj[l]].size(); i++ )
          vout << upf_r[i] << " " << 0.5*vps[iproj[l]][i] << endl;
        vout << endl << endl;
      }
    }
    vout << endl << endl;
    vout.close();
#endif

    // Generate QSO file

    // output potential in QSO format
    cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    cout << "<fpmd:species xmlns:fpmd=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"" << endl;
    cout << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
    cout << "  xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0"  << endl;
    cout << "  species.xsd\">" << endl;
    cout << "<description>" << endl;
    cout << "Translated from UPF format by upf2qso " << release
         << " on " << isodate() << endl;
    cout << upf_pp_info;
    cout << "</description>" << endl;
    cout << "<symbol>" << upf_symbol << "</symbol>" << endl;
    cout << "<atomic_number>" << atomic_number << "</atomic_number>" << endl;
    cout << "<mass>" << mass << "</mass>" << endl;
    cout << "<norm_conserving_pseudopotential>" << endl;
    cout << "<valence_charge>" << upf_zval << "</valence_charge>" << endl;
    cout << "<lmax>" << qso_lmax << "</lmax>" << endl;
    cout << "<llocal>" << upf_llocal << "</llocal>" << endl;
    cout << "<nquad>0</nquad>" << endl;
    cout << "<rquad>0.0</rquad>" << endl;
    cout << "<mesh_spacing>" << mesh_spacing << "</mesh_spacing>" << endl;

    cout.setf(ios::scientific,ios::floatfield);
    if ( upf_nlcc_flag == "T" )
    {
      cout << "<core_density size=\"" << nplin << "\">" << endl;
      for ( int i = 0; i < nplin; i++ )
        cout << setprecision(10) << nlcc_lin[i] << endl;
      cout << "</core_density>" << endl;
    }

    for ( int l = 0; l <= qso_lmax; l++ )
    {
      cout << "<projector l=\"" << l << "\" size=\"" << nplin << "\">"
           << endl;
      cout << "<radial_potential>" << endl;
      if ( iproj[l] == -1 )
      {
        // l == llocal
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(10) << vloc_lin[i] << endl;
      }
      else
      {
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(10) << vps_lin[iproj[l]][i] << endl;
      }
      cout << "</radial_potential>" << endl;
      // find index j corresponding to angular momentum l
      int j = 0;
      while ( upf_wf_l[j] != l && j < upf_nwf ) j++;
      // check if found
      const bool found = ( j != upf_nwf );
      // print wf only if found
      if ( found )
      {
        cout << "<radial_function>" << endl;
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(10) << wf_lin[j][i] << endl;
        cout << "</radial_function>" << endl;
      }
      cout << "</projector>" << endl;
    }
    cout << "</norm_conserving_pseudopotential>" << endl;
    cout << "</fpmd:species>" << endl;

  }
  else if ( upf_version == 2 )
  {
    // process UPF version 2 potential
    seek_str("<PP_INFO>");
    string upf_pp_info;
    bool done = false;
    while (!done)
    {
      getline(cin,buf);
      is.clear();
      is.str(buf);
      is >> s;
      done = ( s == "</PP_INFO>" );
      if ( !done )
      {
        upf_pp_info += buf + '\n';
      }
    }

    // remove all '<' and '>' characters from the PP_INFO field
    // for XML compatibility
    p = upf_pp_info.find_first_of("<>");
    while ( p != string::npos )
    {
      upf_pp_info[p] = ' ';
      p = upf_pp_info.find_first_of("<>");
    }

    string tag = find_start_element("PP_HEADER");

    // get attribute "element"
    string upf_symbol = get_attr(tag,"element");
    cerr << " upf_symbol: " << upf_symbol << endl;

    // get atomic number and mass
    const int atomic_number = pt.z(upf_symbol);
    const double mass = pt.mass(upf_symbol);

    // check if potential is norm-conserving or semi-local
    string pseudo_type = get_attr(tag,"pseudo_type");
    cerr << " pseudo_type = " << pseudo_type << endl;
    if ( pseudo_type!="NC" && pseudo_type!="SL" )
    {
      cerr << " pseudo_type must be NC or SL" << endl;
      return 1;
    }

    // NLCC flag
    string upf_nlcc_flag = get_attr(tag,"core_correction");
    if ( upf_nlcc_flag == "T" )
    {
      cerr << " Potential includes a non-linear core correction" << endl;
    }
    cerr << " upf_nlcc_flag = " << upf_nlcc_flag << endl;

    // XC functional (add in description)
    string upf_functional = get_attr(tag,"functional");
    // add XC functional information to description
    upf_pp_info += "functional = " + upf_functional + '\n';
    cerr << " upf_functional = " << upf_functional << endl;

    // valence charge
    double upf_zval = 0.0;
    string buf = get_attr(tag,"z_valence");
    is.clear();
    is.str(buf);
    is >> upf_zval;
    cerr << " upf_zval = " << upf_zval << endl;

    // max angular momentum
    int upf_lmax;
    buf = get_attr(tag,"l_max");
    is.clear();
    is.str(buf);
    is >> upf_lmax;
    cerr << " upf_lmax = " << upf_lmax << endl;

    // local angular momentum
    int upf_llocal;
    buf = get_attr(tag,"l_local");
    is.clear();
    is.str(buf);
    is >> upf_llocal;
    cerr << " upf_llocal = " << upf_llocal << endl;

    // number of points in mesh
    int upf_mesh_size;
    buf = get_attr(tag,"mesh_size");
    is.clear();
    is.str(buf);
    is >> upf_mesh_size;
    cerr << " upf_mesh_size = " << upf_mesh_size << endl;

    // number of wavefunctions
    int upf_nwf;
    buf = get_attr(tag,"number_of_wfc");
    is.clear();
    is.str(buf);
    is >> upf_nwf;
    cerr << " upf_nwf = " << upf_nwf << endl;

    // number of projectors
    int upf_nproj;
    buf = get_attr(tag,"number_of_proj");
    is.clear();
    is.str(buf);
    is >> upf_nproj;
    cerr << " upf_nproj = " << upf_nproj << endl;

    vector<int> upf_l(upf_nwf);

    // read mesh
    find_start_element("PP_MESH");
    find_start_element("PP_R");
    vector<double> upf_r(upf_mesh_size);
    for ( int i = 0; i < upf_mesh_size; i++ )
      cin >> upf_r[i];
    find_end_element("PP_R");
    find_start_element("PP_RAB");
    vector<double> upf_rab(upf_mesh_size);
    for ( int i = 0; i < upf_mesh_size; i++ )
      cin >> upf_rab[i];
    find_end_element("PP_RAB");
    find_end_element("PP_MESH");

    find_start_element("PP_LOCAL");
    vector<double> upf_vloc(upf_mesh_size);
    for ( int i = 0; i < upf_mesh_size; i++ )
      cin >> upf_vloc[i];
    find_end_element("PP_LOCAL");

    find_start_element("PP_NONLOCAL");
    vector<vector<double> > upf_vnl;
    upf_vnl.resize(upf_nproj);
    vector<int> upf_proj_l(upf_nproj);

    ostringstream os;
    for ( int j = 0; j < upf_nproj; j++ )
    {
      int index, angular_momentum;
      os.str("");
      os << j+1;
      string element_name = "PP_BETA." + os.str();
      tag = find_start_element(element_name);
      cerr << tag << endl;

      buf = get_attr(tag,"index");
      is.clear();
      is.str(buf);
      is >> index;
      cerr << " index = " << index << endl;

      buf = get_attr(tag,"angular_momentum");
      is.clear();
      is.str(buf);
      is >> angular_momentum;
      cerr << " angular_momentum = " << angular_momentum << endl;

      assert(angular_momentum <= upf_lmax);
      upf_proj_l[index-1] = angular_momentum;

      upf_vnl[j].resize(upf_mesh_size);
      for ( int i = 0; i < upf_mesh_size; i++ )
        cin >> upf_vnl[j][i];

      find_end_element(element_name);
    }

    // compute number of projectors for each l
    // nproj_l[l] is the number of projectors having angular momentum l
    vector<int> nproj_l(upf_lmax+1);
    for ( int l = 0; l <= upf_lmax; l++ )
    {
      nproj_l[l] = 0;
      for ( int ip = 0; ip < upf_nproj; ip++ )
        if ( upf_proj_l[ip] == l ) nproj_l[l]++;
    }

    tag = find_start_element("PP_DIJ");
    int size;
    buf = get_attr(tag,"size");
    is.clear();
    is.str(buf);
    is >> size;
    cerr << "PP_DIJ size = " << size << endl;

    if ( size != upf_nproj*upf_nproj )
    {
      cerr << " Number of non-zero Dij differs from number of projectors"
           << endl;
      return 1;
    }
    int upf_ndij = size;

    vector<double> upf_d(upf_ndij);
    for ( int i = 0; i < upf_ndij; i++ )
    {
      cin >> upf_d[i];
    }
    int imax = sqrt(size+1.e-5);
    assert(imax*imax==size);

    // Check if Dij has non-diagonal elements
    // non-diagonal elements are not supported
    for ( int m = 0; m < imax; m++ )
      for ( int n = 0; n < imax; n++ )
        if ( (m != n) && (upf_d[n*imax+m] != 0.0) )
        {
          cerr << " Non-local Dij has off-diagonal elements" << endl;
          cerr << " m=" << m << " n=" << n << endl;
          return 1;
        }

    find_end_element("PP_DIJ");

    find_end_element("PP_NONLOCAL");

    // make table iproj[l] mapping l to iproj
    // vnl(l) is in vnl[iproj[l]] if iproj[l] > -1
    // vlocal if iproj[llocal] = -1
    vector<int> iproj(upf_lmax+2);
    for ( int l = 0; l <= upf_lmax+1; l++ )
      iproj[l] = -1;
    for ( int j = 0; j < upf_nproj; j++ )
      iproj[upf_proj_l[j]] = j;

    // determine angular momentum of local potential in UPF file
    // upf_llocal is the angular momentum of the local potential
    // increase lmax if there are more projectors than wavefunctions
    int qso_lmax = upf_lmax;
    if (upf_lmax < upf_llocal)
    {
      qso_lmax = upf_lmax+1;
    }

    if ( pseudo_type == "SL" )
    {
      find_start_element("PP_PSWFC");
      vector<vector<double> > upf_wf;
      vector<int> upf_wf_l(upf_nwf);
      upf_wf.resize(upf_nwf);
      for ( int j = 0; j < upf_nwf; j++ )
      {
        int index, l;
        os.str("");
        os << j+1;
        string element_name = "PP_CHI." + os.str();
        tag = find_start_element(element_name);
        cerr << tag << endl;

        buf = get_attr(tag,"index");
        is.clear();
        is.str(buf);
        is >> index;
        cerr << " index = " << index << endl;

        buf = get_attr(tag,"l");
        is.clear();
        is.str(buf);
        is >> l;
        cerr << " l = " << l << endl;

        assert(l <= upf_lmax);
        upf_proj_l[index-1] = l;
        upf_wf[j].resize(upf_mesh_size);
        for ( int i = 0; i < upf_mesh_size; i++ )
          cin >> upf_wf[j][i];
      }
      find_end_element("PP_PSWFC");

      // NLCC
      vector<double> upf_nlcc;
      if ( upf_nlcc_flag == "T" )
      {
        find_start_element("PP_NLCC");
        upf_nlcc.resize(upf_mesh_size);
        for ( int i = 0; i < upf_mesh_size; i++ )
          cin >> upf_nlcc[i];
        find_end_element("PP_NLCC");
      }

      // output original data in file upf.dat
      ofstream upf("upf.dat");
      upf << "# vloc" << endl;
      for ( int i = 0; i < upf_vloc.size(); i++ )
        upf << upf_r[i] << " " << upf_vloc[i] << endl;
      upf << endl << endl;
      for ( int j = 0; j < upf_nproj; j++ )
      {
        upf << "# proj j=" << j << endl;
        for ( int i = 0; i < upf_vnl[j].size(); i++ )
          upf << upf_r[i] << " " << upf_vnl[j][i] << endl;
        upf << endl << endl;
      }
      for ( int j = 0; j < upf_nwf; j++ )
      {
        upf << "# wf j=" << j << endl;
        for ( int i = 0; i < upf_wf[j].size(); i++ )
          upf << upf_r[i] << " " << upf_wf[j][i] << endl;
        upf << endl << endl;
      }
      upf.close();

      // print summary
      cerr << "PP_INFO:" << endl << upf_pp_info << endl;
      cerr << "Element: " << upf_symbol << endl;
       cerr << "NLCC: " << upf_nlcc_flag << endl;
      //cerr << "XC: " << upf_xcf[0] << " " << upf_xcf[1] << " "
      //     << upf_xcf[2] << " " << upf_xcf[3] << endl;
      cerr << "Zv: " << upf_zval << endl;
      cerr << "lmax: " << qso_lmax << endl;
      cerr << "llocal: " << upf_llocal << endl;
      cerr << "nwf: " << upf_nwf << endl;
      cerr << "mesh_size: " << upf_mesh_size << endl;

      // compute delta_vnl[l][i] on the upf log mesh

      // divide the projector function by the wavefunction, except if
      // the wavefunction amplitude is smaller than tol, outside of rcut_divide.
      const double tol = 1.e-5;
      const double rcut_divide = 1.0;
      vector<vector<double> > delta_vnl;
      delta_vnl.resize(upf_nproj);
      for ( int j = 0; j < upf_nproj; j++ )
      {
        delta_vnl[j].resize(upf_wf[j].size());
        for ( int i = 0; i < delta_vnl[j].size(); i++ )
        {
          double num = upf_vnl[j][i];
          double den = upf_wf[upf_proj_l[j]][i];

          delta_vnl[j][i] = 0.0;
          if ( upf_r[i] < rcut_divide )
          {
            // near the origin
            if ( i == 0 && fabs(den) < tol )
            {
              // i = 0 for linear mesh, r = 0.0: use nearest value
              delta_vnl[j][i] = upf_vnl[j][1] / upf_wf[upf_proj_l[j]][1];
            }
            else
            {
              // other points near the origin
              delta_vnl[j][i] = num / den;
            }
          }
          else
          {
            // wavefunction gets small at large r.
            // Assume that delta_vnl is zero when that happens
            if ( fabs(den) > tol )
              delta_vnl[j][i] = num / den;
          }
        }
      }

      vector<vector<double> > vps;
      vps.resize(upf_nproj+1);
      for ( int j = 0; j < upf_nproj; j++ )
      {
        vps[j].resize(upf_mesh_size);
        for ( int i = 0; i < delta_vnl[j].size(); i++ )
          vps[j][i] = upf_vloc[i] + delta_vnl[j][i];
      }

      // interpolate functions on linear mesh
      const double mesh_spacing = 0.01;
      int nplin = (int) (rcut / mesh_spacing);
      vector<double> f(upf_mesh_size), fspl(upf_mesh_size);

      // interpolate NLCC
      vector<double> nlcc_lin(nplin);
      if ( upf_nlcc_flag == "T" )
      {
        assert(upf_mesh_size==upf_nlcc.size());
        for ( int i = 0; i < upf_nlcc.size(); i++ )
          f[i] = upf_nlcc[i];
        int n = upf_nlcc.size();
        int bcnat_left = 0;
        double yp_left = 0.0;
        int bcnat_right = 1;
        double yp_right = 0.0;
        spline(n,&upf_r[0],&f[0],yp_left,yp_right,
               bcnat_left,bcnat_right,&fspl[0]);

        for ( int i = 0; i < nplin; i++ )
        {
          double r = i * mesh_spacing;
          if ( r >= upf_r[0] )
            splint(n,&upf_r[0],&f[0],&fspl[0],r,&nlcc_lin[i]);
          else
            // use value closest to the origin for r=0
            nlcc_lin[i] = upf_nlcc[0];
          if ( fabs(nlcc_lin[i]) < 1.e-12 )
            nlcc_lin[i] = 0.0;
        }
      }

      // interpolate vloc
      // factor 0.5: convert from Ry in UPF to Hartree in QSO
      for ( int i = 0; i < upf_vloc.size(); i++ )
        f[i] = 0.5 * upf_vloc[i];

      int n = upf_vloc.size();
      int bcnat_left = 0;
      double yp_left = 0.0;
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);

      vector<double> vloc_lin(nplin);
      for ( int i = 0; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        if ( r >= upf_r[0] )
          splint(n,&upf_r[0],&f[0],&fspl[0],r,&vloc_lin[i]);
        else
          // use value closest to the origin for r=0
          vloc_lin[i] = 0.5 * upf_vloc[0];
      }

      // interpolate vps[j], j=0, nproj-1
      vector<vector<double> > vps_lin;
      vps_lin.resize(vps.size());
      for ( int j = 0; j < vps.size(); j++ )
      {
        vps_lin[j].resize(nplin);
      }

      for ( int j = 0; j < upf_nproj; j++ )
      {
        // factor 0.5: convert from Ry in UPF to Hartree in QSO
        for ( int i = 0; i < upf_vloc.size(); i++ )
          f[i] = 0.5 * vps[j][i];

        int n = upf_vloc.size();
        int bcnat_left = 0;
        double yp_left = 0.0;
        int bcnat_right = 1;
        double yp_right = 0.0;
        spline(n,&upf_r[0],&f[0],yp_left,yp_right,
               bcnat_left,bcnat_right,&fspl[0]);

        for ( int i = 0; i < nplin; i++ )
        {
          double r = i * mesh_spacing;
          if ( r >= upf_r[0] )
            splint(n,&upf_r[0],&f[0],&fspl[0],r,&vps_lin[j][i]);
          else
            vps_lin[j][i] = 0.5 * vps[j][0];
        }
      }

      // write potentials in gnuplot format on file vlin.dat
      ofstream vlin("vlin.dat");
      for ( int l = 0; l <= qso_lmax; l++ )
      {
        vlin << "# v, l=" << l << endl;
        if ( iproj[l] == -1 )
        {
          // l == llocal
          for ( int i = 0; i < nplin; i++ )
            vlin << i*mesh_spacing << " " << vloc_lin[i] << endl;
          vlin << endl << endl;
        }
        else
        {
          for ( int i = 0; i < nplin; i++ )
            vlin << i*mesh_spacing << " " << vps_lin[iproj[l]][i] << endl;
          vlin << endl << endl;
        }
      }

      // interpolate wavefunctions on the linear mesh

      vector<vector<double> > wf_lin;
      wf_lin.resize(upf_nwf);
      for ( int j = 0; j < upf_nwf; j++ )
      {
        wf_lin[j].resize(nplin);
        assert(upf_wf[j].size()<=f.size());
        for ( int i = 0; i < upf_wf[j].size(); i++ )
        {
          if ( upf_r[i] > 0.0 )
            f[i] = upf_wf[j][i] / upf_r[i];
          else
          {
            // value at origin, depending on angular momentum
            if ( upf_wf_l[j] == 0 )
            {
              // l=0: take value closest to r=0
              f[i] = upf_wf[j][1]/upf_r[1];
            }
            else
            {
              // l>0:
              f[i] = 0.0;
            }
          }
        }

        int n = upf_wf[j].size();
        // choose boundary condition at origin depending on angular momentum
        int bcnat_left = 0;
        double yp_left = 0.0;
        if ( upf_wf_l[j] == 1 )
        {
          bcnat_left = 1; // use natural bc
          double yp_left = 0.0; // not used
        }
        int bcnat_right = 1;
        double yp_right = 0.0;
        spline(n,&upf_r[0],&f[0],yp_left,yp_right,
               bcnat_left,bcnat_right,&fspl[0]);

        for ( int i = 0; i < nplin; i++ )
        {
          double r = i * mesh_spacing;
          if ( r >= upf_r[0] )
            splint(n,&upf_r[0],&f[0],&fspl[0],r,&wf_lin[j][i]);
          else
          {
            // r < upf_r[0]
            assert(upf_r[0]>0.0);
            // compute value near origin, depending on angular momentum
            if ( upf_wf_l[j] == 0 )
            {
              // l=0: take value closest to r=0
              wf_lin[j][i] = upf_wf[j][0]/upf_r[0];
            }
            else
            {
              // l>0:
              wf_lin[j][i] = upf_wf[j][0] * r / ( upf_r[0] * upf_r[0] );
            }
          }
        }

        vlin << "# phi, l=" << upf_l[j] << endl;
        for ( int i = 0; i < nplin; i++ )
          vlin << i*mesh_spacing << " " << wf_lin[j][i] << endl;
        vlin << endl << endl;
      }

      cerr << " interpolation done" << endl;

  #if 1
      // output potential on log mesh
      ofstream vout("v.dat");
      for ( int l = 0; l <= qso_lmax; l++ )
      {
        vout << "# v, l=" << l << endl;
        if ( iproj[l] == -1 )
        {
          // l == llocal
          for ( int i = 0; i < upf_vloc.size(); i++ )
            vout << upf_r[i] << " " << 0.5*upf_vloc[i] << endl;
          vout << endl << endl;
        }
        else
        {
          for ( int i = 0; i < vps[iproj[l]].size(); i++ )
            vout << upf_r[i] << " " << 0.5*vps[iproj[l]][i] << endl;
          vout << endl << endl;
        }
      }
      vout << endl << endl;
      vout.close();
  #endif

      // Generate QSO file

      // output potential in QSO format
      cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
      cout << "<fpmd:species xmlns:fpmd=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"" << endl;
      cout << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
      cout << "  xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0"  << endl;
      cout << "  species.xsd\">" << endl;
      cout << "<description>" << endl;
      cout << "Translated from UPF format by upf2qso " << release
           << " on " << isodate() << endl;
      cout << upf_pp_info;
      cout << "</description>" << endl;
      cout << "<symbol>" << upf_symbol << "</symbol>" << endl;
      cout << "<atomic_number>" << atomic_number << "</atomic_number>" << endl;
      cout << "<mass>" << mass << "</mass>" << endl;
      cout << "<norm_conserving_pseudopotential>" << endl;
      cout << "<valence_charge>" << upf_zval << "</valence_charge>" << endl;
      cout << "<lmax>" << qso_lmax << "</lmax>" << endl;
      cout << "<llocal>" << upf_llocal << "</llocal>" << endl;
      cout << "<nquad>0</nquad>" << endl;
      cout << "<rquad>0.0</rquad>" << endl;
      cout << "<mesh_spacing>" << mesh_spacing << "</mesh_spacing>" << endl;

      cout.setf(ios::scientific,ios::floatfield);
      if ( upf_nlcc_flag == "T" )
      {
        cout << "<core_density size=\"" << nplin << "\">" << endl;
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(10) << nlcc_lin[i] << endl;
        cout << "</core_density>" << endl;
      }

      for ( int l = 0; l <= qso_lmax; l++ )
      {
        cout << "<projector l=\"" << l << "\" size=\"" << nplin << "\">"
             << endl;
        cout << "<radial_potential>" << endl;
        if ( iproj[l] == -1 )
        {
          // l == llocal
          for ( int i = 0; i < nplin; i++ )
            cout << setprecision(10) << vloc_lin[i] << endl;
        }
        else
        {
          for ( int i = 0; i < nplin; i++ )
            cout << setprecision(10) << vps_lin[iproj[l]][i] << endl;
        }
        cout << "</radial_potential>" << endl;
        // find index j corresponding to angular momentum l
        int j = 0;
        while ( upf_wf_l[j] != l && j < upf_nwf ) j++;
        // check if found
        const bool found = ( j != upf_nwf );
        // print wf only if found
        if ( found )
        {
          cout << "<radial_function>" << endl;
          for ( int i = 0; i < nplin; i++ )
            cout << setprecision(10) << wf_lin[j][i] << endl;
          cout << "</radial_function>" << endl;
        }
        cout << "</projector>" << endl;
      }
      cout << "</norm_conserving_pseudopotential>" << endl;
      cout << "</fpmd:species>" << endl;
    } // if SL

    if ( pseudo_type == "NC" )
    {
      // NLCC
      vector<double> upf_nlcc;
      if ( upf_nlcc_flag == "T" )
      {
        find_start_element("PP_NLCC");
        upf_nlcc.resize(upf_mesh_size);
        for ( int i = 0; i < upf_mesh_size; i++ )
          cin >> upf_nlcc[i];
        find_end_element("PP_NLCC");
      }

      cerr << " NC potential" << endl;
      // output original data in file upf.dat
      ofstream upf("upf.dat");
      upf << "# vloc" << endl;
      for ( int i = 0; i < upf_vloc.size(); i++ )
        upf << upf_r[i] << " " << upf_vloc[i] << endl;
      upf << endl << endl;
      for ( int j = 0; j < upf_nproj; j++ )
      {
        upf << "# proj j=" << j << endl;
        for ( int i = 0; i < upf_vnl[j].size(); i++ )
          upf << upf_r[i] << " " << upf_vnl[j][i] << endl;
        upf << endl << endl;
      }

      upf << "# dij " << endl;
      for ( int j = 0; j < upf_d.size(); j++ )
      {
        upf << j << " " << upf_d[j] << endl;
      }
      upf.close();

      // print summary
      cerr << "PP_INFO:" << endl << upf_pp_info << endl;
      cerr << "Element: " << upf_symbol << endl;
      cerr << "NLCC: " << upf_nlcc_flag << endl;
      // cerr << "XC: " << upf_xcf[0] << " " << upf_xcf[1] << " "
      //      << upf_xcf[2] << " " << upf_xcf[3] << endl;
      cerr << "Zv: " << upf_zval << endl;
      cerr << "lmax: " << qso_lmax << endl;
      cerr << "nproj: " << upf_nproj << endl;
      cerr << "mesh_size: " << upf_mesh_size << endl;

      // interpolate functions on linear mesh
      const double mesh_spacing = 0.01;
      int nplin = (int) (rcut / mesh_spacing);
      vector<double> f(upf_mesh_size), fspl(upf_mesh_size);

      // interpolate NLCC
      vector<double> nlcc_lin(nplin);
      if ( upf_nlcc_flag == "T" )
      {
        assert(upf_mesh_size==upf_nlcc.size());
        for ( int i = 0; i < upf_nlcc.size(); i++ )
          f[i] = upf_nlcc[i];
        int n = upf_nlcc.size();
        int bcnat_left = 0;
        double yp_left = 0.0;
        int bcnat_right = 1;
        double yp_right = 0.0;
        spline(n,&upf_r[0],&f[0],yp_left,yp_right,
               bcnat_left,bcnat_right,&fspl[0]);

        for ( int i = 0; i < nplin; i++ )
        {
          double r = i * mesh_spacing;
          if ( r >= upf_r[0] )
            splint(n,&upf_r[0],&f[0],&fspl[0],r,&nlcc_lin[i]);
          else
            // use value closest to the origin for r=0
            nlcc_lin[i] = upf_nlcc[0];
          if ( fabs(nlcc_lin[i]) < 1.e-12 )
            nlcc_lin[i] = 0.0;
        }
      }

      // interpolate vloc
      // factor 0.5: convert from Ry in UPF to Hartree in QSO
      for ( int i = 0; i < upf_vloc.size(); i++ )
        f[i] = 0.5 * upf_vloc[i];

      int n = upf_vloc.size();
      int bcnat_left = 0;
      double yp_left = 0.0;
      int bcnat_right = 1;
      double yp_right = 0.0;
      spline(n,&upf_r[0],&f[0],yp_left,yp_right,
             bcnat_left,bcnat_right,&fspl[0]);

      vector<double> vloc_lin(nplin);
      for ( int i = 0; i < nplin; i++ )
      {
        double r = i * mesh_spacing;
        if ( r >= upf_r[0] )
          splint(n,&upf_r[0],&f[0],&fspl[0],r,&vloc_lin[i]);
        else
          // use value closest to the origin for r=0
          vloc_lin[i] = 0.5 * upf_vloc[0];
      }

      // interpolate vnl[j], j=0, nproj-1
      vector<vector<double> > vnl_lin;
      vnl_lin.resize(upf_nproj);
      for ( int j = 0; j < vnl_lin.size(); j++ )
      {
        vnl_lin[j].resize(nplin);
      }

      for ( int j = 0; j < upf_nproj; j++ )
      {
        // interpolate projectors
        // note: upf_vnl contains r*projector
        // See UPF documentation at http://www.quantum-espresso.org/
        //   pseudopotentials/unified-pseudopotential-format
        assert(f.size()>=upf_vnl[j].size());
        for ( int i = 0; i < upf_vnl[j].size(); i++ )
          f[i] = upf_vnl[j][i];

        int n = f.size();
        int bcnat_left = 1;
        double yp_left = 0.0;
        int bcnat_right = 1;
        double yp_right = 0.0;
        spline(n,&upf_r[0],&f[0],yp_left,yp_right,
               bcnat_left,bcnat_right,&fspl[0]);

        for ( int i = 0; i < nplin; i++ )
        {
          double r = i * mesh_spacing;
          if ( r >= upf_r[0] )
            splint(n,&upf_r[0],&f[0],&fspl[0],r,&vnl_lin[j][i]);
          else
            vnl_lin[j][i] = upf_vnl[j][0];
          if ( fabs(vnl_lin[j][i]) < 1.e-12 )
            vnl_lin[j][i] = 0.0;
        }
      }

      // write local potential and projectors in gnuplot format on file vlin.dat
      ofstream vlin("vlin.dat");
      vlin << "# vlocal" << endl;
      for ( int i = 0; i < nplin; i++ )
        vlin << vloc_lin[i] << endl;
      vlin << endl << endl;
      for ( int iproj = 0; iproj < vnl_lin.size(); iproj++ )
      {
        vlin << "# projector, l=" << upf_proj_l[iproj] << endl;
        for ( int i = 0; i < nplin; i++ )
          vlin << i*mesh_spacing << " " << vnl_lin[iproj][i] << endl;
        vlin << endl << endl;
      }

      // Generate QSO file

      // output potential in QSO format
      cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
      cout << "<fpmd:species xmlns:fpmd=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"" << endl;
      cout << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
      cout << "  xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0"  << endl;
      cout << "  species.xsd\">" << endl;
      cout << "<description>" << endl;
      cout << "Translated from UPF format by upf2qso " << release
           << " on " << isodate() << endl;
      cout << upf_pp_info;
      cout << "</description>" << endl;
      cout << "<symbol>" << upf_symbol << "</symbol>" << endl;
      cout << "<atomic_number>" << atomic_number << "</atomic_number>" << endl;
      cout << "<mass>" << mass << "</mass>" << endl;
      cout << "<norm_conserving_semilocal_pseudopotential>" << endl;
      cout << "<valence_charge>" << upf_zval << "</valence_charge>" << endl;
      cout << "<mesh_spacing>" << mesh_spacing << "</mesh_spacing>" << endl;

      cout.setf(ios::scientific,ios::floatfield);
      if ( upf_nlcc_flag == "T" )
      {
        cout << "<core_density size=\"" << nplin << "\">" << endl;
        for ( int i = 0; i < nplin; i++ )
          cout << setprecision(10) << nlcc_lin[i] << endl;
        cout << "</core_density>" << endl;
      }

      // local potential
      cout << "<local_potential size=\"" << nplin << "\">" << endl;
      for ( int i = 0; i < nplin; i++ )
        cout << setprecision(10) << vloc_lin[i] << endl;
      cout << "</local_potential>" << endl;

      // projectors
      // note: vnl_lin contains r[i]*projector[i]
      int ip = 0;
      for ( int l = 0; l <= upf_lmax; l++ )
      {
        for ( int i = 0; i < nproj_l[l]; i++ )
        {
          cout << "<projector l=\"" << l << "\" i=\""
               << i+1 << "\" size=\"" << nplin << "\">"
               << endl;
          // value at r=0:
          // quadratic extrapolation of vnl_lin(r)/r to r=0 if l==0
          if ( l == 0 )
          {
            // use quadratic extrapolation to zero
            const double h = mesh_spacing;
            const double v = (4.0/3.0)*vnl_lin[ip][1]/h -
                             (1.0/3.0)*vnl_lin[ip][2]/(2*h);
            cout << setprecision(10) << v << endl;
          }
          else
          {
            cout << setprecision(10) << 0.0 << endl;
          }
          for ( int j = 1; j < nplin; j++ )
          {
            const double r = j * mesh_spacing;
            cout << setprecision(10) << vnl_lin[ip][j]/r << endl;
          }
          ip++;
          cout << "</projector>" << endl;
        }
      }

      // d_ij
      int ibase = 0;
      int jbase = 0;
      for ( int l = 0; l <= upf_lmax; l++ )
      {
        for ( int i = 0; i < upf_nproj; i++ )
          for ( int j = 0; j < upf_nproj; j++ )
          {
            if ( (upf_proj_l[i] == l) && (upf_proj_l[j] == l) )
            {
              int ij = i + j*upf_nproj;
              cout << "<d_ij l=\"" << l << "\""
                   << " i=\"" << i-ibase+1 << "\" j=\"" << j-jbase+1
                   << "\"> " << setprecision(10) << 0.5*upf_d[ij] << " </d_ij>"
                   << endl;
            }
          }
        ibase += nproj_l[l];
        jbase += nproj_l[l];
      }

      cout << "</norm_conserving_semilocal_pseudopotential>" << endl;
      cout << "</fpmd:species>" << endl;
    }
  } // version 1 or 2
  return 0;
}
