////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
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
//  KPGridConnectivity.C
//
////////////////////////////////////////////////////////////////////////////////

#include "KPGridConnectivity.h"
using namespace std;

#define TAG_Overlaps_first_kx 6
#define TAG_Overlaps_first_ky 7
#define TAG_Overlaps_first_kz 8
#define TAG_Overlaps_second_kx 9
#define TAG_Overlaps_second_ky 10
#define TAG_Overlaps_second_kz 11
#define TAG_Overlaps_local 12

////////////////////////////////////////////////////////////////////////////////
KPConnectivity::KPConnectivity(const Sample& s)
{
  const double kdist_tol = 1.e-5;
  int onpe0 = ( ( s.ctxt_.myrow()==0 ) && ( s.ctxt_.mycol()==0 ) );

#ifdef DEBUG
  if ( onpe0 )
    cout << "Building kp grid connection for quad HF exchange convergence\n";
#endif

  // get number of kpoints
  nkpoints_ = s.wf.nkp();
  comm_     = s.wf.sd(0,0)->context().comm();

  // get maximum number of local states
  int nStateLoc_= s.wf.sd(0,0)->nstloc();
  MPI_Allreduce(&nStateLoc_,&nStateMax_,1,MPI_INT,MPI_MAX,comm_);

  // allocate memory for kpoint weights
  weight_.resize(nkpoints_);

  // allocate memory for the list of integrals
  integral_kx_.resize(nkpoints_);
  integral_ky_.resize(nkpoints_);
  integral_kz_.resize(nkpoints_);

  // allocate memory for the list of neighbours
  first_neighbour_kx_.resize(nkpoints_);
  first_neighbour_ky_.resize(nkpoints_);
  first_neighbour_kz_.resize(nkpoints_);
  second_neighbour_kx_.resize(nkpoints_);
  second_neighbour_ky_.resize(nkpoints_);
  second_neighbour_kz_.resize(nkpoints_);

  // allocate memory for the symmetry tags
  first_symmetric_kx_.resize(nkpoints_);
  first_symmetric_ky_.resize(nkpoints_);
  first_symmetric_kz_.resize(nkpoints_);
  second_symmetric_kx_.resize(nkpoints_);
  second_symmetric_ky_.resize(nkpoints_);
  second_symmetric_kz_.resize(nkpoints_);

  // allocate memory for the index of overlaps
  local_ig_overlap_.resize(nkpoints_);
  first_ig_overlap_kx_.resize(nkpoints_);
  first_ig_overlap_ky_.resize(nkpoints_);
  first_ig_overlap_kz_.resize(nkpoints_);
  second_ig_overlap_kx_.resize(nkpoints_);
  second_ig_overlap_ky_.resize(nkpoints_);
  second_ig_overlap_kz_.resize(nkpoints_);

  // allocate memory for the index of overlaps
  first_T_overlap_kx_.resize(nkpoints_);
  first_T_overlap_ky_.resize(nkpoints_);
  first_T_overlap_kz_.resize(nkpoints_);
  second_T_overlap_kx_.resize(nkpoints_);
  second_T_overlap_ky_.resize(nkpoints_);
  second_T_overlap_kz_.resize(nkpoints_);

  // allocate memory for the distances
  first_distance_kx_.resize(nkpoints_);
  first_distance_ky_.resize(nkpoints_);
  first_distance_kz_.resize(nkpoints_);
  second_distance_kx_.resize(nkpoints_);
  second_distance_ky_.resize(nkpoints_);
  second_distance_kz_.resize(nkpoints_);

  // allocate memory for the overlaps
  local_overlap_.resize(nkpoints_*nStateMax_);
  first_overlap_kx_.resize(nkpoints_*nStateMax_);
  first_overlap_ky_.resize(nkpoints_*nStateMax_);
  first_overlap_kz_.resize(nkpoints_*nStateMax_);
  second_overlap_kx_.resize(nkpoints_*nStateMax_);
  second_overlap_ky_.resize(nkpoints_*nStateMax_);
  second_overlap_kz_.resize(nkpoints_*nStateMax_);

  // allocate memory for the communications
  local_send_buff_.resize(nStateMax_);
  first_send_buff_kx_.resize(nStateMax_);
  first_send_buff_ky_.resize(nStateMax_);
  first_send_buff_kz_.resize(nStateMax_);
  second_send_buff_kx_.resize(nStateMax_);
  second_send_buff_ky_.resize(nStateMax_);
  second_send_buff_kz_.resize(nStateMax_);

  // get the dimensionality of the kpoint grid
  DimX_=0;
  DimY_=0;
  DimZ_=0;
  for ( int iKp=0 ; iKp<nkpoints_ ; iKp++ )
  {
    // test direct difference
    D3vector dk =   s.wf.kpoint(iKp) - s.wf.kpoint(0);
    if ( abs(dk.x)>kdist_tol ) DimX_=1;
    if ( abs(dk.y)>kdist_tol ) DimY_=1;
    if ( abs(dk.z)>kdist_tol ) DimZ_=1;

    // test diference with symmetric
    dk =   s.wf.kpoint(iKp) + s.wf.kpoint(0);
    if ( abs(dk.x)>kdist_tol ) DimX_=1;
    if ( abs(dk.y)>kdist_tol ) DimY_=1;
    if ( abs(dk.z)>kdist_tol ) DimZ_=1;
  }

  // get the length of each vector of the
  // reciprocal cell
  double length_kx=length(s.wf.cell().b(0));
  double length_ky=length(s.wf.cell().b(1));
  double length_kz=length(s.wf.cell().b(2));

  // get the connectivity of the grid
  for ( int iKpi=0 ; iKpi<nkpoints_ ; iKpi++ )
  {
    // initialize min distance
    double dminx1=1.e10;
    double dminy1=1.e10;
    double dminz1=1.e10;
    double dminx2=1.e10;
    double dminy2=1.e10;
    double dminz2=1.e10;

    // initialize indices to -1
    first_neighbour_kx_[iKpi]=-1;
    first_neighbour_ky_[iKpi]=-1;
    first_neighbour_kz_[iKpi]=-1;
    second_neighbour_kx_[iKpi]=-1;
    second_neighbour_ky_[iKpi]=-1;
    second_neighbour_kz_[iKpi]=-1;

    // initialize distances
    first_distance_kx_[iKpi]=1;
    first_distance_ky_[iKpi]=1;
    first_distance_kz_[iKpi]=1;
    second_distance_kx_[iKpi]=1;
    second_distance_ky_[iKpi]=1;
    second_distance_kz_[iKpi]=1;

    // if this kpoint is part of the integration grid
    if ( s.wf.weight(iKpi)!=0.0 )
    {
      // explore the remaining k points
      ConnectivityComplete_=1;
      for ( int iKpj=0 ; iKpj<nkpoints_ ; iKpj++ )
      {
        // if this kpoint is part of the integration grid
        if ( s.wf.weight(iKpj)!=0.0 )
        {
          for ( int itx = -1; itx <=1; itx++ )
          for ( int ity = -1; ity <=1; ity++ )
          for ( int itz = -1; itz <=1; itz++ )
          {
	    D3vector T(itx, ity, itz);
            // first test direct difference
            D3vector dk = s.wf.kpoint(iKpi) - s.wf.kpoint(iKpj) - T;

            // make sure that we are not just considering 0 vector
            if ( dk.x!=0 || dk.y!=0 || dk.z!=0 )
            {
              // kx direction
              if ( abs(dk.x)<dminx1+kdist_tol && abs(dk.y)<kdist_tol &&
                   abs(dk.z)<kdist_tol )
              {
                // copy first neighbour in second
                dminx2=dminx1;
                second_neighbour_kx_[iKpi]=first_neighbour_kx_[iKpi];
                second_symmetric_kx_[iKpi]=first_symmetric_kx_[iKpi];
                second_T_overlap_kx_[iKpi]=first_T_overlap_kx_[iKpi];
                second_distance_kx_[iKpi] =first_distance_kx_[iKpi];
                // set new first neighbour
                dminx1=abs(dk.x);
                first_neighbour_kx_[iKpi]=iKpj;
                first_symmetric_kx_[iKpi]=0;
                first_T_overlap_kx_[iKpi]=T;
                first_distance_kx_[iKpi] =dk.x*length_kx;
              }
              else if ( abs(dk.x)<dminx2+kdist_tol && abs(dk.y)<kdist_tol &&
                        abs(dk.z)<kdist_tol )
              {
                // set new second neighbour
                dminx2=abs(dk.x);
                second_neighbour_kx_[iKpi]=iKpj;
                second_symmetric_kx_[iKpi]=0;
                second_T_overlap_kx_[iKpi]=T;
                second_distance_kx_[iKpi] =dk.x*length_kx;
              }

              // ky direction
              if ( abs(dk.x)<kdist_tol && abs(dk.y)<dminy1+kdist_tol &&
                   abs(dk.z)<kdist_tol )
              {
                // copy first neighbour in second
                dminy2=dminy1;
                second_neighbour_ky_[iKpi]=first_neighbour_ky_[iKpi];
                second_symmetric_ky_[iKpi]=first_symmetric_ky_[iKpi];
                second_T_overlap_ky_[iKpi]=first_T_overlap_ky_[iKpi];
                second_distance_ky_[iKpi] =first_distance_ky_[iKpi];
                // set new first neighbour
                dminy1=abs(dk.y);
                first_neighbour_ky_[iKpi]=iKpj;
                first_symmetric_ky_[iKpi]=0;
                first_T_overlap_ky_[iKpi]=T;
                first_distance_ky_[iKpi] =dk.y*length_ky;
              }
              else if ( abs(dk.x)<kdist_tol && abs(dk.y)<dminy2+kdist_tol &&
                        abs(dk.z)<kdist_tol )
              {
                // set new second neighbour
                dminy2=abs(dk.y);
                second_neighbour_ky_[iKpi]=iKpj;
                second_symmetric_ky_[iKpi]=0;
                second_T_overlap_ky_[iKpi]=T;
                second_distance_ky_[iKpi] =dk.y*length_ky;
              }

              // kz direction
              if ( abs(dk.x)<kdist_tol && abs(dk.y)<kdist_tol &&
                   abs(dk.z)<dminz1+kdist_tol )
              {
                // copy first neighbour in second
                dminz2=dminz1;
                second_neighbour_kz_[iKpi]=first_neighbour_kz_[iKpi];
                second_symmetric_kz_[iKpi]=first_symmetric_kz_[iKpi];
                second_T_overlap_kz_[iKpi]=first_T_overlap_kz_[iKpi];
                second_distance_kz_[iKpi] =first_distance_kz_[iKpi];
                // set new first neighbour
                dminz1=abs(dk.z);
                first_neighbour_kz_[iKpi]=iKpj;
                first_symmetric_kz_[iKpi]=0;
                first_T_overlap_kz_[iKpi]=T;
                first_distance_kz_[iKpi] =dk.z*length_kz;
              }
              else if ( abs(dk.x)<kdist_tol && abs(dk.y)<kdist_tol &&
                        abs(dk.z)<dminz2+kdist_tol )
              {
                // set new second neighbour
                dminz2=abs(dk.z);
                second_neighbour_kz_[iKpi]=iKpj;
                second_symmetric_kz_[iKpi]=0;
                second_T_overlap_kz_[iKpi]=T;
                second_distance_kz_[iKpi] =dk.z*length_kz;
              }
            }

            // then test difference with symmetric
            // (except for gamma and 0.5 vectors)
            if ( ( s.wf.kpoint(iKpj).x!=0.0 &&
                   fabs(s.wf.kpoint(iKpj).x)!=0.5 ) ||
                 ( s.wf.kpoint(iKpj).y!=0.0 &&
                   fabs(s.wf.kpoint(iKpj).y)!=0.5 ) ||
                 ( s.wf.kpoint(iKpj).z!=0.0 &&
                   fabs(s.wf.kpoint(iKpj).z)!=0.5 ) )
            {
              dk = s.wf.kpoint(iKpi) + s.wf.kpoint(iKpj) + T;

              // make sure that we are not just considering 0 vector
              if ( dk.x!=0 || dk.y!=0 || dk.z!=0 )
              {
                // kx direction
                if ( abs(dk.x)<dminx1+kdist_tol && abs(dk.y)<kdist_tol &&
                     abs(dk.z)<kdist_tol )
                {
                  // copy first neighbour in second
                  dminx2=dminx1;
                  second_neighbour_kx_[iKpi]=first_neighbour_kx_[iKpi];
                  second_symmetric_kx_[iKpi]=first_symmetric_kx_[iKpi];
                  second_T_overlap_kx_[iKpi]=first_T_overlap_kx_[iKpi];
                  second_distance_kx_[iKpi] =first_distance_kx_[iKpi];
                  // set new first neighbour
                  dminx1=abs(dk.x);
                  first_neighbour_kx_[iKpi]=iKpj;
                  first_symmetric_kx_[iKpi]=1;
                  first_T_overlap_kx_[iKpi]=T;
                  first_distance_kx_[iKpi] =dk.x*length_kx;
                }
                else if ( abs(dk.x)<dminx2+kdist_tol && abs(dk.y)<kdist_tol &&
                          abs(dk.z)<kdist_tol )
                {
                  // set new second neighbour
                  dminx2=abs(dk.x);
                  second_neighbour_kx_[iKpi]=iKpj;
                  second_symmetric_kx_[iKpi]=1;
                  second_T_overlap_kx_[iKpi]=T;
                  second_distance_kx_[iKpi] =dk.x*length_kx;
                }

                // ky direction
                if ( abs(dk.x)<kdist_tol && abs(dk.y)<dminy1+kdist_tol &&
                     abs(dk.z)<kdist_tol )
                {
                  // copy first neighbour in second
                  dminy2=dminy1;
                  second_neighbour_ky_[iKpi]=first_neighbour_ky_[iKpi];
                  second_symmetric_ky_[iKpi]=first_symmetric_ky_[iKpi];
                  second_T_overlap_ky_[iKpi]=first_T_overlap_ky_[iKpi];
                  second_distance_ky_[iKpi] =first_distance_ky_[iKpi];
                  // set new first neighbour
                  dminy1=abs(dk.y);
                  first_neighbour_ky_[iKpi]=iKpj;
                  first_symmetric_ky_[iKpi]=1;
                  first_T_overlap_ky_[iKpi]=T;
                  first_distance_ky_[iKpi] =dk.y*length_ky;
                }
                else if ( abs(dk.x)<kdist_tol && abs(dk.y)<dminy2+kdist_tol &&
                          abs(dk.z)<kdist_tol )
                {
                  // set new second neighbour
                  dminy2=abs(dk.y);
                  second_neighbour_ky_[iKpi]=iKpj;
                  second_symmetric_ky_[iKpi]=1;
                  second_T_overlap_ky_[iKpi]=T;
                  second_distance_ky_[iKpi] =dk.y*length_ky;
                }

                // kz direction
                if ( abs(dk.x)<kdist_tol && abs(dk.y)<kdist_tol &&
                     abs(dk.z)<dminz1+kdist_tol )
                {
                  // copy first neighbour in second
                  dminz2=dminz1;
                  second_neighbour_kz_[iKpi]=first_neighbour_kz_[iKpi];
                  second_symmetric_kz_[iKpi]=first_symmetric_kz_[iKpi];
                  second_T_overlap_kz_[iKpi]=first_T_overlap_kz_[iKpi];
                  second_distance_kz_[iKpi] =first_distance_kz_[iKpi];
                  // set new first neighbour
                  dminz1=abs(dk.z);
                  first_neighbour_kz_[iKpi]=iKpj;
                  first_symmetric_kz_[iKpi]=1;
                  first_T_overlap_kz_[iKpi]=T;
                  first_distance_kz_[iKpi] =dk.z*length_kz;
                }
                else if ( abs(dk.x)<kdist_tol && abs(dk.y)<kdist_tol &&
                          abs(dk.z)<dminz2+kdist_tol )
                {
                  // set new second neighbour
                  dminz2=abs(dk.z);
                  second_neighbour_kz_[iKpi]=iKpj;
                  second_symmetric_kz_[iKpi]=1;
                  second_T_overlap_kz_[iKpi]=T;
                  second_distance_kz_[iKpi] =dk.z*length_kz;
                }
              }
            }
          }
        }
      }

      // check if we found the connectivity for this kpoint
      if ( first_neighbour_kx_[iKpi]==-1 || second_neighbour_kx_[iKpi]==-1 )
      {
        if ( onpe0 )
          cout << "\nWarning: cannot find connectivity kx for kpoint "
               << iKpi << "\n\n";
        ConnectivityComplete_=0;
        return;
      }
      if ( first_neighbour_ky_[iKpi]==-1 || second_neighbour_ky_[iKpi]==-1 )
      {
        if ( onpe0 )
          cout << "\nWarning: cannot find connectivity ky for kpoint "
               << iKpi << "\n\n";
        ConnectivityComplete_=0;
        return;
      }
      if ( first_neighbour_kz_[iKpi]==-1 || second_neighbour_kz_[iKpi]==-1 )
      {
        if ( onpe0 )
          cout << "\nWarning: cannot find connectivity kz for kpoint "
               << iKpi << "\n\n";
        ConnectivityComplete_=0;
        return;
      }

      // check if the grid is irregular (not implemented yet)
      // take into account only the kpoint with non zero weight
      if (fabs(second_distance_kx_[iKpi]+first_distance_kx_[iKpi])>kdist_tol ||
          fabs(second_distance_ky_[iKpi]+first_distance_ky_[iKpi])>kdist_tol ||
          fabs(second_distance_kz_[iKpi]+first_distance_kz_[iKpi])>kdist_tol)
      {
        if ( onpe0 )
        {
          cout << " Found an irregular grid at kpoint "
               << iKpi << ": cannot use connectivity.\n";
          cout << " distances kx: " << second_distance_kx_[iKpi]
               << " " << first_distance_kx_[iKpi] << " neigbour kx:"
               << second_neighbour_kx_[iKpi] << " "
               << first_neighbour_kx_[iKpi] << "\n";
          cout << " distances ky: " << second_distance_ky_[iKpi]
               << " " << first_distance_ky_[iKpi] << " neigbour ky:"
               << second_neighbour_ky_[iKpi] << " "
               << first_neighbour_ky_[iKpi] << "\n";
          cout << " distances kz: " << second_distance_kz_[iKpi]
               << " " << first_distance_kz_[iKpi] << " neigbour kz:"
               << second_neighbour_kz_[iKpi] << " "
               << first_neighbour_kz_[iKpi] << "\n";
        }
        ConnectivityComplete_=0;
        return;
      }
    }

    // compute the weight of this kpoint
    // this weight should be non zero only for kpoints
    // taken into account during the computation
    if ( s.wf.weight(iKpi)!=0.0 )
    {
      weight_[iKpi]=fabs((first_distance_kx_[iKpi]-second_distance_kx_[iKpi]) *
                         (first_distance_ky_[iKpi]-second_distance_ky_[iKpi]) *
                         (first_distance_ky_[iKpi]-second_distance_ky_[iKpi]))/
                        (8.0 * length_kx * length_ky * length_kz);
    }
    else
    {
      weight_[iKpi]=0.0;
    }

#ifdef DEBUG
    if ( onpe0 )
      cout << "Weight kpoint " << iKpi << " = " <<  weight_[iKpi] << "\n";
#endif
  }

  // compute the value of the integral of v^2/(x^2+y^2+z^2) for v = {x,y,z},
  // normed to a unitary volume of integration (deltaX*deltaY*deltaZ = 1).
  double deltaX=fabs(first_distance_kx_[0]-second_distance_kx_[0]);
  double deltaY=fabs(first_distance_ky_[0]-second_distance_ky_[0]);
  double deltaZ=fabs(first_distance_kz_[0]-second_distance_kz_[0]);
  double volume=deltaX*deltaY*deltaZ;
  deltaX=deltaX/pow(volume,3.0);
  deltaY=deltaY/pow(volume,3.0);
  deltaZ=deltaZ/pow(volume,3.0);

  double xmin=-deltaX/2.0;
  double xmax= deltaX/2.0;
  double ymin=-deltaY/2.0;
  double ymax= deltaY/2.0;
  double zmin=-deltaZ/2.0;
  double zmax= deltaZ/2.0;

  // estimate numerically the integral
  // start by creating a grid
  int n=1000;
  vector<double> x(n);
  vector<double> y(n);
  vector<double> z(n);

  for ( int i=0 ; i<n ; i++ )
  {
    x[i]=xmin+deltaX/(n+1)*i;
    y[i]=ymin+deltaY/(n+1)*i;
    z[i]=zmin+deltaZ/(n+1)*i;
  }

  // integrate analytically over tha first coordinate
  // then integrate numerically
  double IntX=0.0;
  double IntY=0.0;
  double IntZ=0.0;

  for ( int i=0 ; i<n ; i++ )
  for ( int j=0 ; j<n ; j++ )
  {
    double SQRTx2py2=sqrt(x[i]*x[i]+y[j]*y[j]);
    double SQRTy2pz2=sqrt(y[i]*y[i]+z[j]*z[j]);
    double SQRTx2pz2=sqrt(x[i]*x[i]+z[j]*z[j]);
    double xinf=xmin/SQRTy2pz2;
    double xsup=xmax/SQRTy2pz2;
    double yinf=ymin/SQRTx2pz2;
    double ysup=ymax/SQRTx2pz2;
    double zinf=zmin/SQRTx2py2;
    double zsup=zmax/SQRTx2py2;
    IntX+=SQRTy2pz2*(atan(xsup)-atan(xinf));
    IntY+=SQRTx2pz2*(atan(ysup)-atan(yinf));
    IntZ+=SQRTx2py2*(atan(zsup)-atan(zinf));
  }
  IntX/=(n*n);
  IntY/=(n*n);
  IntZ/=(n*n);

  // the sum of the three integral should be equal to 2
  // renormalize
  //!! note: this step should not be necessary
  //!! replace with an assert to check if the sum is correct
  double sum = IntX + IntY + IntZ;
  IntX *= ( 2.0 / sum );
  IntY *= ( 2.0 / sum );
  IntZ *= ( 2.0 / sum );

  // the final value of the integrals is  1 - IntV
  IntX = 1.0 - IntX;
  IntY = 1.0 - IntY;
  IntZ = 1.0 - IntZ;

#ifdef DEBUG
  if ( onpe0 )
    cout << "Curvature integrals: "
         << IntX << " , " << IntY << " , " << IntZ << "\n";
#endif

  // we associate these values with all kpoints
  // (based on the assumption of a regular grid)
  for ( int iKpi=0 ; iKpi<nkpoints_ ; iKpi++ )
  {
    integral_kx_[iKpi]=IntX;
    integral_ky_[iKpi]=IntY;
    integral_kz_[iKpi]=IntZ;
  }

#ifdef DEBUG
  if ( onpe0 )
    cout << "Succesfully built kpoint grid connection: the kpoint sampling is "
         << DimX_+DimY_+DimZ_ << "D\n";
#endif
}

////////////////////////////////////////////////////////////////////////////////
KPConnectivity::~KPConnectivity()
{}

////////////////////////////////////////////////////////////////////////////////
void KPConnectivity::SetOverlapIndices(Basis *vbasis_)
{
  // first, find the indices corresponding to the
  // first translations of the reciprocal cell
  int indices[3][3][3];
  // init the indices to -1
  for ( int i=-1 ; i<=1 ; i++ )
  for ( int j=-1 ; j<=1 ; j++ )
  for ( int k=-1 ; k<=1 ; k++ )
  {
    indices[i+1][j+1][k+1]=-1;
  }
  // find the corresponding indices in the G basis
  const int ngloc = vbasis_->localsize();
  const UnitCell& cell = vbasis_->cell();
  for ( int ig=0 ; ig < ngloc; ig++ )
  {
    // get the corresponding index
    //!! D3vector G;
    //!! G.x=vbasis_->gx(ig+G_local_basis_size*0);
    //!! G.y=vbasis_->gx(ig+G_local_basis_size*1);
    //!! G.z=vbasis_->gx(ig+G_local_basis_size*2);
    //!! int ix = (int) round( G * cell.a(0) / 6.283185307179586 );
    //!! int iy = (int) round( G * cell.a(1) / 6.283185307179586 );
    //!! int iz = (int) round( G * cell.a(2) / 6.283185307179586 );
    int ix = vbasis_->idx(3*ig+0);
    int iy = vbasis_->idx(3*ig+1);
    int iz = vbasis_->idx(3*ig+2);
    if ( -1 <= ix && ix <= 1 &&
         -1 <= iy && iy <= 1 &&
         -1 <= iz && iz <= 1 )
      indices[1+ix][1+iy][1+iz]=ig;
  }
  // set the indices of the overlap for each kpoint
  for ( int iKpi=0 ; iKpi<nkpoints_ ; iKpi++ )
  {
    // local overlap
    local_ig_overlap_[iKpi] = indices[1][1][1];

    // first neighbors
    if ( !( iKpi>first_neighbour_kx_[iKpi] || first_symmetric_kx_[iKpi] ) )
    {
      int ix = -(int)first_T_overlap_kx_[iKpi].x;
      int iy = -(int)first_T_overlap_kx_[iKpi].y;
      int iz = -(int)first_T_overlap_kx_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      first_ig_overlap_kx_[iKpi] = indices[ix+1][iy+1][iz+1];
    }
    else
    {
      int ix =  (int)first_T_overlap_kx_[iKpi].x;
      int iy =  (int)first_T_overlap_kx_[iKpi].y;
      int iz =  (int)first_T_overlap_kx_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      first_ig_overlap_kx_[iKpi] = indices[ix+1][iy+1][iz+1];
    }

    if ( !( iKpi>first_neighbour_ky_[iKpi] || first_symmetric_ky_[iKpi] ) )
    {
      int ix = -(int)first_T_overlap_ky_[iKpi].x;
      int iy = -(int)first_T_overlap_ky_[iKpi].y;
      int iz = -(int)first_T_overlap_ky_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      first_ig_overlap_ky_[iKpi] = indices[ix+1][iy+1][iz+1];
    }
    else
    {
      int ix =  (int)first_T_overlap_ky_[iKpi].x;
      int iy =  (int)first_T_overlap_ky_[iKpi].y;
      int iz =  (int)first_T_overlap_ky_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      first_ig_overlap_ky_[iKpi] = indices[ix+1][iy+1][iz+1];
    }

    if ( !( iKpi>first_neighbour_kz_[iKpi] || first_symmetric_kz_[iKpi] ) )
    {
      int ix = -(int)first_T_overlap_kz_[iKpi].x;
      int iy = -(int)first_T_overlap_kz_[iKpi].y;
      int iz = -(int)first_T_overlap_kz_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      first_ig_overlap_kz_[iKpi] = indices[ix+1][iy+1][iz+1];
    }
    else
    {
      int ix =  (int)first_T_overlap_kz_[iKpi].x;
      int iy =  (int)first_T_overlap_kz_[iKpi].y;
      int iz =  (int)first_T_overlap_kz_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      first_ig_overlap_kz_[iKpi] = indices[ix+1][iy+1][iz+1];
    }

    // second neighbors
    if ( !( iKpi>second_neighbour_kx_[iKpi] || second_symmetric_kx_[iKpi] ) )
    {
      int ix = -(int)second_T_overlap_kx_[iKpi].x;
      int iy = -(int)second_T_overlap_kx_[iKpi].y;
      int iz = -(int)second_T_overlap_kx_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      second_ig_overlap_kx_[iKpi] = indices[ix+1][iy+1][iz+1];
    }
    else
    {
      int ix =  (int)second_T_overlap_kx_[iKpi].x;
      int iy =  (int)second_T_overlap_kx_[iKpi].y;
      int iz =  (int)second_T_overlap_kx_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      second_ig_overlap_kx_[iKpi] = indices[ix+1][iy+1][iz+1];
    }

    if ( !( iKpi>second_neighbour_ky_[iKpi] || second_symmetric_ky_[iKpi] ) )
    {
      int ix = -(int)second_T_overlap_ky_[iKpi].x;
      int iy = -(int)second_T_overlap_ky_[iKpi].y;
      int iz = -(int)second_T_overlap_ky_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      second_ig_overlap_ky_[iKpi] = indices[ix+1][iy+1][iz+1];
    }
    else
    {
      int ix =  (int)second_T_overlap_ky_[iKpi].x;
      int iy =  (int)second_T_overlap_ky_[iKpi].y;
      int iz =  (int)second_T_overlap_ky_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      second_ig_overlap_ky_[iKpi] = indices[ix+1][iy+1][iz+1];
    }

    if ( !( iKpi>second_neighbour_kz_[iKpi] || second_symmetric_kz_[iKpi] ) )
    {
      int ix = -(int)second_T_overlap_kz_[iKpi].x;
      int iy = -(int)second_T_overlap_kz_[iKpi].y;
      int iz = -(int)second_T_overlap_kz_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      second_ig_overlap_kz_[iKpi] = indices[ix+1][iy+1][iz+1];
    }
    else
    {
      int ix =  (int)second_T_overlap_kz_[iKpi].x;
      int iy =  (int)second_T_overlap_kz_[iKpi].y;
      int iz =  (int)second_T_overlap_kz_[iKpi].z;
      assert(ix*ix<2 && iy*iy<2 && iz*iz<2);
      second_ig_overlap_kz_[iKpi] = indices[ix+1][iy+1][iz+1];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void KPConnectivity::InitOverlaps()
{
  for ( int iKpi=0 ; iKpi<nkpoints_ ; iKpi++ )
  {
    for ( int iState=0 ; iState<nStateMax_ ; iState++ )
    {
      local_overlap_[iKpi*nStateMax_+iState]=0.0;
      first_overlap_kx_[iKpi*nStateMax_+iState]=0.0;
      first_overlap_ky_[iKpi*nStateMax_+iState]=0.0;
      first_overlap_kz_[iKpi*nStateMax_+iState]=0.0;
      second_overlap_kx_[iKpi*nStateMax_+iState]=0.0;
      second_overlap_ky_[iKpi*nStateMax_+iState]=0.0;
      second_overlap_kz_[iKpi*nStateMax_+iState]=0.0;
    }
  }
}

// add overlaps:
// accumulate the sum of the square modulus of the
// overlaps over the first neighbours
// take occupation number into account
////////////////////////////////////////////////////////////////////////////////
void KPConnectivity::AddOverlap(int iKpi, int iKpj, int iLocStatei,
  complex<double> *valueDirect, complex<double> *valueSymmetric,
  double occupation)
{
  // fast return
  if ( occupation==0.0 ) return;

  // local overlap
  if ( iKpi==iKpj && local_ig_overlap_[iKpi] >= 0 )
    local_overlap_[iKpi*nStateMax_+iLocStatei] +=
      norm(valueDirect[ local_ig_overlap_[iKpi] ])  * occupation;

  // case of kpx neighbour
  if ( first_neighbour_kx_[iKpi]==iKpj && first_ig_overlap_kx_[iKpi]>=0 )
  {
    // store overlap value for both states
    first_overlap_kx_[iKpi*nStateMax_+iLocStatei] +=
      ( first_symmetric_kx_[iKpi] ?
        norm(valueSymmetric[ first_ig_overlap_kx_[iKpi] ]) :
        norm(valueDirect[ first_ig_overlap_kx_[iKpi] ]) ) * occupation;
  }

  if ( second_neighbour_kx_[iKpi]==iKpj && second_ig_overlap_kx_[iKpi]>=0 )
  {
    // store overlap value for both states
    second_overlap_kx_[iKpi*nStateMax_+iLocStatei] +=
      ( second_symmetric_kx_[iKpi] ?
        norm(valueSymmetric[ second_ig_overlap_kx_[iKpi] ]) :
        norm(valueDirect[ second_ig_overlap_kx_[iKpi] ]) ) * occupation;
  }

  // case of kpy neighbour
  if ( first_neighbour_ky_[iKpi]==iKpj && first_ig_overlap_ky_[iKpi]>=0 )
  {
    // store overlap value for both states
    first_overlap_ky_[iKpi*nStateMax_+iLocStatei] +=
      ( first_symmetric_ky_[iKpi] ?
        norm(valueSymmetric[ first_ig_overlap_ky_[iKpi] ]) :
        norm(valueDirect[ first_ig_overlap_ky_[iKpi] ]) ) * occupation;
  }

  if ( second_neighbour_ky_[iKpi]==iKpj && second_ig_overlap_ky_[iKpi]>=0 )
  {
    // store overlap value for both states
    second_overlap_ky_[iKpi*nStateMax_+iLocStatei] +=
      ( second_symmetric_ky_[iKpi] ?
        norm(valueSymmetric[ second_ig_overlap_ky_[iKpi] ]) :
        norm(valueDirect[ second_ig_overlap_ky_[iKpi] ]) ) * occupation;
  }

  // case of kpz neighbour
  if ( first_neighbour_kz_[iKpi]==iKpj && first_ig_overlap_kz_[iKpi]>=0 )
  {
    // store overlap value for both states
    first_overlap_kz_[iKpi*nStateMax_+iLocStatei] +=
      ( first_symmetric_kz_[iKpi] ?
        norm(valueSymmetric[ first_ig_overlap_kz_[iKpi] ]) :
        norm(valueDirect[ first_ig_overlap_kz_[iKpi] ]) ) * occupation;
  }

  if ( second_neighbour_kz_[iKpi]==iKpj && second_ig_overlap_kz_[iKpi]>=0 )
  {
    // store overlap value for both states
    second_overlap_kz_[iKpi*nStateMax_+iLocStatei] +=
      ( second_symmetric_kz_[iKpi] ?
        norm(valueSymmetric[ second_ig_overlap_kz_[iKpi] ]) :
        norm(valueDirect[ second_ig_overlap_kz_[iKpi] ]) ) * occupation;
  }
}
#undef norm

// array permutation
////////////////////////////////////////////////////////////////////////////////
void KPConnectivity::StartPermutation(int iKp, int iSendTo, int iRecvFr)
{
  // local overlaps
  if (1)
  {
    // copy overlaps into send buffer
    for ( int i=0 ; i<nStateMax_ ; i++ )
    {
      local_send_buff_[i]=local_overlap_[iKp*nStateMax_+i];
    }
    // send the data in non blocking mode
    MPI_Isend((void *) &local_send_buff_[0], nStateMax_, MPI_DOUBLE,
      iSendTo, TAG_Overlaps_local, comm_, &send_request_local_ );
    // receive the data in non blocking mode
    MPI_Irecv((void *) &local_overlap_[iKp*nStateMax_], nStateMax_, MPI_DOUBLE,
      iRecvFr, TAG_Overlaps_local, comm_, &recv_request_local_ );
  }

  // if direction kx is used
  if (DimX_==1)
  {
    // copy overlaps into send buffer
    for ( int i=0 ; i<nStateMax_ ; i++ )
    {
      first_send_buff_kx_[i]=first_overlap_kx_[iKp*nStateMax_+i];
      second_send_buff_kx_[i]=second_overlap_kx_[iKp*nStateMax_+i];
    }
    // send the data in non blocking mode
    MPI_Isend((void *) &first_send_buff_kx_[0], nStateMax_, MPI_DOUBLE,
      iSendTo, TAG_Overlaps_first_kx, comm_, &send_request_first_kx_ );
    MPI_Isend((void *) &second_send_buff_kx_[0], nStateMax_, MPI_DOUBLE,
      iSendTo, TAG_Overlaps_second_kx, comm_, &send_request_second_kx_ );
    // receive the data in non blocking mode
    MPI_Irecv((void *) &first_overlap_kx_[iKp*nStateMax_], nStateMax_,
      MPI_DOUBLE, iRecvFr, TAG_Overlaps_first_kx, comm_,
      &recv_request_first_kx_ );
    MPI_Irecv((void *) &second_overlap_kx_[iKp*nStateMax_], nStateMax_,
      MPI_DOUBLE, iRecvFr, TAG_Overlaps_second_kx, comm_,
      &recv_request_second_kx_ );
  }

  // if direction ky is used
  if (DimY_==1)
  {
    // copy overlaps into send buffer
    for ( int i=0 ; i<nStateMax_ ; i++ )
    {
      first_send_buff_ky_[i]=first_overlap_ky_[iKp*nStateMax_+i];
      second_send_buff_ky_[i]=second_overlap_ky_[iKp*nStateMax_+i];
    }
    // send the data in non blocking mode
    MPI_Isend((void *) &first_send_buff_ky_[0], nStateMax_, MPI_DOUBLE,
      iSendTo, TAG_Overlaps_first_ky, comm_, &send_request_first_ky_ );
    MPI_Isend((void *) &second_send_buff_ky_[0], nStateMax_, MPI_DOUBLE,
      iSendTo, TAG_Overlaps_second_ky, comm_, &send_request_second_ky_ );
    // receive the data in non blocking mode
    MPI_Irecv((void *) &first_overlap_ky_[iKp*nStateMax_], nStateMax_,
      MPI_DOUBLE, iRecvFr, TAG_Overlaps_first_ky, comm_,
      &recv_request_first_ky_ );
    MPI_Irecv((void *) &second_overlap_ky_[iKp*nStateMax_], nStateMax_,
      MPI_DOUBLE, iRecvFr, TAG_Overlaps_second_ky, comm_,
      &recv_request_second_ky_ );
  }

  // if direction kz is used
  if (DimZ_==1)
  {
    // copy overlaps into send buffer
    for ( int i=0 ; i<nStateMax_ ; i++ )
    {
      first_send_buff_kz_[i]=first_overlap_kz_[iKp*nStateMax_+i];
      second_send_buff_kz_[i]=second_overlap_kz_[iKp*nStateMax_+i];
    }
    // send the data in non blocking mode
    MPI_Isend((void *) &first_send_buff_kz_[0], nStateMax_, MPI_DOUBLE,
      iSendTo, TAG_Overlaps_first_kz, comm_, &send_request_first_kz_ );
    MPI_Isend((void *) &second_send_buff_kz_[0], nStateMax_, MPI_DOUBLE,
      iSendTo, TAG_Overlaps_second_kz, comm_, &send_request_second_kz_ );
    // receive the data in non blocking mode
    MPI_Irecv((void *) &first_overlap_kz_[iKp*nStateMax_], nStateMax_,
      MPI_DOUBLE, iRecvFr, TAG_Overlaps_first_kz, comm_,
      &recv_request_first_kz_ );
    MPI_Irecv((void *) &second_overlap_kz_[iKp*nStateMax_], nStateMax_,
      MPI_DOUBLE, iRecvFr, TAG_Overlaps_second_kz, comm_,
      &recv_request_second_kz_ );
  }
}

////////////////////////////////////////////////////////////////////////////////
void KPConnectivity::EndPermutation()
{
  // local overlaps
  if (1)
  {
    MPI_Status status_recv_local;
    MPI_Status status_send_local;
    MPI_Wait( &recv_request_local_, &status_recv_local);
    MPI_Wait( &send_request_local_, &status_send_local);
  }

  // if direction kx is used
  if (DimX_==1)
  {
    MPI_Status status_recv_first;
    MPI_Status status_send_first;
    MPI_Status status_recv_second;
    MPI_Status status_send_second;
    MPI_Wait( &recv_request_first_kx_, &status_recv_first);
    MPI_Wait( &send_request_first_kx_, &status_send_first);
    MPI_Wait( &recv_request_second_kx_, &status_recv_second);
    MPI_Wait( &send_request_second_kx_, &status_send_second);
  }

  // if direction ky is used
  if (DimY_==1)
  {
    MPI_Status status_recv_first;
    MPI_Status status_send_first;
    MPI_Status status_recv_second;
    MPI_Status status_send_second;
    MPI_Wait( &recv_request_first_ky_, &status_recv_first);
    MPI_Wait( &send_request_first_ky_, &status_send_first);
    MPI_Wait( &recv_request_second_ky_, &status_recv_second);
    MPI_Wait( &send_request_second_ky_, &status_send_second);
  }

  // if direction kz is used
  if (DimZ_==1)
  {
    MPI_Status status_recv_first;
    MPI_Status status_send_first;
    MPI_Status status_recv_second;
    MPI_Status status_send_second;
    MPI_Wait( &recv_request_first_kz_, &status_recv_first);
    MPI_Wait( &send_request_first_kz_, &status_send_first);
    MPI_Wait( &recv_request_second_kz_, &status_recv_second);
    MPI_Wait( &send_request_second_kz_, &status_send_second);
  }
}

// access to the accumulated overlaps
////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::overlaps_local(int iKp, int iLocState)
{
  return local_overlap_[iKp*nStateMax_+iLocState];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::overlaps_first_kx(int iKp, int iLocState)
{
  return first_overlap_kx_[iKp*nStateMax_+iLocState];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::overlaps_second_kx(int iKp, int iLocState)
{
  return second_overlap_kx_[iKp*nStateMax_+iLocState];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::overlaps_first_ky(int iKp, int iLocState)
{
  return first_overlap_ky_[iKp*nStateMax_+iLocState];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::overlaps_second_ky(int iKp, int iLocState)
{
  return second_overlap_ky_[iKp*nStateMax_+iLocState];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::overlaps_first_kz(int iKp, int iLocState)
{
  return first_overlap_kz_[iKp*nStateMax_+iLocState];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::overlaps_second_kz(int iKp, int iLocState)
{
  return second_overlap_kz_[iKp*nStateMax_+iLocState];
}

// access to the neighbour distances
////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::distance_first_kx(int iKp)
{
  return first_distance_kx_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::distance_second_kx(int iKp)
{
  return second_distance_kx_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::distance_first_ky(int iKp)
{
  return first_distance_ky_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::distance_second_ky(int iKp)
{
  return second_distance_ky_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::distance_first_kz(int iKp)
{
  return first_distance_kz_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::distance_second_kz(int iKp)
{
  return second_distance_kz_[iKp];
}

// access to the integral values
////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::integral_kx(int iKp)
{
  return integral_kx_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::integral_ky(int iKp)
{
  return integral_ky_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::integral_kz(int iKp)
{
  return integral_kz_[iKp];
}

////////////////////////////////////////////////////////////////////////////////
double KPConnectivity::weight(int iKp)
{
  return weight_[iKp];
}
