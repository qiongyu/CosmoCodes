#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <fftw3.h>
#include <complex>
#include <omp.h>

using namespace std;

int Ngrid = 1024;
long Ngrid3 = Ngrid*Ngrid*Ngrid;
double Lbox = 8000;
int nbin = 80;
double completeness = 1.;

struct _gal
{
  double x, y, z, weight, nz;
};

double exp_x(double x)
{
  double y;

  if(fabs(x) < 0.001)
  {
    double x2 = x*x;
    y = 1. + x + x2/2. + x2*x/6. + x2*x2/24.;
  }
  else y = exp(x);

  return y;
}

double sinc(double x)
{
  if(fabs(x) < 1e-8) return 1.;
  else return sin(x)/x;
}

template <class T>
T* initialize_vector(long dim)
{

  T *x;
  x = new T[dim];

  for(long i=0;i<dim;i++) x[i] = 0;

  return x;
}

template <class V>
V*** initialize_3D_array(long dim1, long dim2, long dim3)
{

  V ***x;

  x = new V**[dim1];

  for(long i=0;i<dim1;i++)
  {
    x[i] = new V*[dim2];
    for(long j=0;j<dim2;j++)
    {
      x[i][j] = new V[dim3];
      for(long k=0;k<dim3;k++) x[i][j][k] = 0.0;
    }
  }

  return x;
}

double* initialize_fft_grid_vectors(double dk)
{

  double *x;
  x = new double[Ngrid];

  for(int i=0;i<Ngrid;i++)
  {
    if(i < Ngrid/2) x[i] = i*dk;
    else x[i] = (i-Ngrid)*dk;
  }

  return x;
}

double* initialize_fft_grid_vectors2()
{

  double *x;
  x = new double[Ngrid];

  for(int i=0;i<Ngrid;i++)
  {
    if(i < Ngrid/2) x[i] = i*i;
    else x[i] = (Ngrid-i)*(Ngrid-i);
  }

  return x;
}

template <class U>
int kill_3D_array(U ***x, long dim1, long dim2)
{
  for(long i=0;i<dim1;i++)
  {
    for(long j=0;j<dim2;j++)
    {
      delete[] x[i][j];
    }
    delete[] x[i];
  }
  delete[] x;

  return 0;
}

template <class K>
K calc_mean(K ***x)
{
  K sum=0;
  for(long i=0;i<Ngrid;i++)
  {
    for(long j=0;j<Ngrid;j++)
    {
      for(long k=0;k<Ngrid;k++) sum += x[i][j][k];
    }
  }
  sum = sum/Ngrid3;
  return sum;
}

long get_num_lines(string infile)
{
  long nline=0;
  string line;
  ifstream ff;

  ff.open(infile.c_str(), ios::in);
  if(ff.is_open())
  {
    while(!ff.eof())
    {
      getline(ff,line);
      nline++;
    }
    ff.close();
  }

//It always counts one extra line for some reason 
//(only EOF indicator on last line??)
  return nline-1;
}

int cic_mass(_gal *p, long N, double ***rho, double box)
{
  double lside = Ngrid/box;
  double lside3 = pow(Ngrid/box,3);

  for(long i=0;i<N;i++)
  {
//Because the box is periodic, the Ngrid+1th grid point is the 0th grid point!
//That's why Ngrid == Ncell!
    double dd1 = p[i].x*lside;
    double dd2 = p[i].y*lside;
    double dd3 = p[i].z*lside;
    long i1 = floor(dd1);
    long i2 = floor(dd2);
    long i3 = floor(dd3);
    long j1 = i1 + 1;
    long j2 = i2 + 1;
    long j3 = i3 + 1;
    dd1 = dd1 - i1;
    dd2 = dd2 - i2;
    dd3 = dd3 - i3;
    double de1 = 1. - dd1;
    double de2 = 1. - dd2;
    double de3 = 1. - dd3;
//If the upper bound of the cell is the Ngrid+1th grid point, then set it back
//to the 0th grid point!
    if(j1 == Ngrid) j1 = 0;
    if(j2 == Ngrid) j2 = 0;
    if(j3 == Ngrid) j3 = 0;
    assert(i1 >= 0 && i1 < Ngrid &&
           i2 >= 0 && i2 < Ngrid &&
           i3 >= 0 && i3 < Ngrid &&
           j1 >= 0 && j1 < Ngrid &&
           j2 >= 0 && j2 < Ngrid &&
           j3 >= 0 && j3 < Ngrid);
/*    if(i1 < 0 || i1 > Ngrid ||
           i2 < 0 || i2 > Ngrid || 
           i3 < 0 || i3 > Ngrid ||
           j1 < 0 || j1 > Ngrid || 
           j2 < 0 || j2 > Ngrid ||
           j3 < 0 || j3 > Ngrid){
    cout << i1 << ' ' << i2 << ' ' << i3 << endl;
    cout << j1 << ' ' << j2 << ' ' << j3 << endl;
    cout << i << ' ' << p[i].x << ' ' << p[i].y << ' ' << p[i].z << endl; }*/

//Need to divide by (Lbox/Ngrid)^3 to get a number density (otherwise you
//are just counting the _TOTAL_ number of galaxies in the cell!)
    rho[i1][i2][i3] += de1*de2*de3 *p[i].weight *lside3;
    rho[j1][i2][i3] += dd1*de2*de3 *p[i].weight *lside3;
    rho[i1][j2][i3] += de1*dd2*de3 *p[i].weight *lside3;
    rho[j1][j2][i3] += dd1*dd2*de3 *p[i].weight *lside3;
    rho[i1][i2][j3] += de1*de2*dd3 *p[i].weight *lside3;
    rho[j1][i2][j3] += dd1*de2*dd3 *p[i].weight *lside3;
    rho[i1][j2][j3] += de1*dd2*dd3 *p[i].weight *lside3;
    rho[j1][j2][j3] += dd1*dd2*dd3 *p[i].weight *lside3;
  }
/*  cout << rho[0][0][0] << endl;
  cout << rho[Ncell][Ncell][Ncell] << endl;
  exit(0);*/

  return 0;
}

int binpk(double ***pk3d, double *kgrid2, double *km, double *pk,
          double norm, double shot)
{
//Bin in k...
  cout << "Done with FFT...on to binning using " << nbin << " bins..." << endl;
  double bin_width = log10(Ngrid/2.)/nbin;

  double *size_bin, *tot_lk, *tot_kpk;
  size_bin = initialize_vector<double>(nbin);
  tot_lk = initialize_vector<double>(nbin);
  tot_kpk = initialize_vector<double>(nbin); //stores k*P(k)

//Because our bins only go up to Ngrid/2 we only need to consider
//P(i,j,k) where 0 < i,j,k < Ngrid/2. 
//However, we do need to consider 0 > -i,-j,-k > -Ngrid/2.
//Luckily, P(-i,-j,-k) = P(Ngrid-i, Ngrid-j, Ngrid-k) which corresponds
//to P(i,j,k) for grid/2 <= i,j,k < Ngrid...
//So we can loop in i,j,k from 0 to Ngrid.
  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        if(i==0 && j==0 && k==0) continue;

        double d2 = kgrid2[i] + kgrid2[j] + kgrid2[k];
//The above is just kx*kx + ky*ky + kz*kz
//Corresponds to:
//i*i+j*j+k*k IF i,j,k<Ngrid/2
//(Ngrid-i)*(Ngrid-i)+(Ngrid-j)*(Ngrid-j)*(Ngrid-k)*(Ngrid-k) IF i,j,k>Ngrid/2
//where the second case represents the -i,-j,-k points of interest.
        double ld = 0.5*log10(d2);
        int bin = floor(ld/bin_width);

        if(bin >= 0 && bin < nbin)
        {
          size_bin[bin]++;
          tot_lk[bin] += ld;
          tot_kpk[bin] += (sqrt(d2) * pk3d[i][j][k]);
        }
      }
    }
  }

  for(int i=0;i<nbin;i++)
  {
    if(size_bin[i] > 0)
    {
      km[i] = pow(10.,tot_lk[i]/size_bin[i]);
      pk[i] = tot_kpk[i]/size_bin[i] / km[i] / norm - shot;

      km[i] = km[i] * (2.*M_PI/Lbox);
//      pk[i] = pk[i] * pow(Lbox,3);
//      cout << i << ' ' << size_bin[i] << endl;
    }
  }

  delete[] size_bin;
  delete[] tot_lk;
  delete[] tot_kpk;

  return 0;
}

int powspec(double ***delta, int nbin, double box, double norm, double shot,
            double *km, double *pk)
{
//Set up grid values in k...(don't need this for fft but need for binning)
  double *kgrid2;
  kgrid2 = initialize_fft_grid_vectors2();
  double dk = 2.*M_PI/Lbox;
  double *kgrid;
  kgrid = initialize_fft_grid_vectors(dk);

//Do fft to get P(k)...
  double *delta_in;
  long r2csize=Ngrid*Ngrid*(Ngrid+2);
  delta_in = initialize_vector<double>(r2csize);
  int count=0;

  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        delta_in[count] = delta[i][j][k];
        count++;
//For r2c, array padding (2 elements) at the end of each k-row. So skip 2 if
//k = Ngrid-1 (last element of k-row).
      }
      count += 2;
    }
  }
  assert(count == r2csize);

  cout << "Calculating the power spectrum (FFT-ing)..." << endl;

  fftw_init_threads();
  fftw_plan_with_nthreads(4);

  fftw_plan plan;
//Do fft in-place!
  plan = fftw_plan_dft_r2c_3d(Ngrid,Ngrid,Ngrid,delta_in,
                             (fftw_complex *)delta_in, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup_threads();

  double lcico2;
  lcico2 = 0.5*Lbox/Ngrid;

  double ***pk_3D;
  pk_3D = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  
  long r2csizek = Ngrid/2+1;
  double lside3 = pow(box/Ngrid,3);
  #pragma omp parallel for schedule(dynamic) 
  for(int i=0;i<Ngrid;i++)
  {
    double sinckx,sincky,sinckz,wcic;
    sinckx = sinc(kgrid[i]*lcico2);
    sinckx *= sinckx;
    for(int j=0;j<Ngrid;j++)
    {
      sincky = sinc(kgrid[j]*lcico2);
      sincky *= sincky;
      for(int k=0;k<r2csizek;k++)
      {
        long dex = (i*Ngrid + j)*r2csizek + k;
        long redex = 2*dex;
        long imdex = redex+1;
        sinckz = sinc(kgrid[k]*lcico2);
        sinckz *= sinckz;
        wcic = sinckx*sincky*sinckz;
        //for some reason I need to divide by Ngrid**3...this appears to be due
        //to a difference between how the FFT is defined in fftw versus IDL...
        //not sure which one is correct????
        //delta_in[redex] = delta_in[redex]/Ngrid3;
        //delta_in[imdex] = delta_in[imdex]/Ngrid3;
        delta_in[redex] = delta_in[redex]*lside3 /wcic;
        delta_in[imdex] = delta_in[imdex]*lside3 /wcic;
        //hopefully this multiplies by the complex conjugate...
        pk_3D[i][j][k] = delta_in[redex]*delta_in[redex]+
                         + delta_in[imdex]*delta_in[imdex];
//For k=(1,...,Ngrid/2-1) we have:
//delta[ks=(Ngrid-k)] = conjugate(delta[k])
//So:
//P[i][j][ks] = delta[i][j][ks]*conjugate(delta[i][j][ks])
//            = delta[i][j][Ngrid-k]*conjugate(delta[i][j][Ngrid-ks])
//            = conjugate(delta[k])*conjugate(conjugate(delta[k]))
//            = conjugate(delta[k])*delta[k]
//            = P[i][j][k]
//k=0 and k=Ngrid/2 are special cases 
//(but note that P[i][j][Ngrid-Ngrid/2] = P[i][j][Ngrid/2], hence, below we
//only need to treat the k=0 case separately)
        if(k!=0) pk_3D[i][j][Ngrid-k] = pk_3D[i][j][k];
      }
    }
  }
  delete[] delta_in;
  delete[] kgrid;

  binpk(pk_3D,kgrid2,km,pk,norm,shot);
  
  delete[] kgrid2;
  kill_3D_array(pk_3D,Ngrid,Ngrid);

  return 0;
}

_gal* read_file(string fname, long *N)
{
  //Read in galaxies
  ifstream infile;

  cout << "Reading " << fname << endl;
  *N = get_num_lines(fname);
  cout << "Number of objects = " << *N << endl;

  _gal *g;
  g = new _gal[*N];
  infile.open( fname.c_str(), ios::in );

  for(long i=0;i<*N;i++)
  {
    string line;
    double temp;

    getline(infile,line);
    stringstream ss(line);
    ss >> g[i].x >> g[i].y >> g[i].z >> g[i].weight >> g[i].nz;
//    g[i].nz *= completeness;
//    g[i].weight = 1;
  }
  cout << g[0].x << ' ' << g[0].y << ' ' << g[0].z << ' ' 
       << g[0].weight << endl;

  return g;
}

int main(int argc, char* argv[])
{
  int ii=1;
  string in, randin;
  while(ii < argc)
  {
    string arg = argv[ii];
    if(arg == "-d")
    {
      ii++;
      in = argv[ii];
      ii++;
    }
    if(arg == "-r")
    {
      ii++;
      randin = argv[ii];
      ii++;
    }
  }
      
  string ifile = in+".txt";
  string randifile = randin+".txt";
  cout << "Using " << Ngrid << " grid points..." << endl;

  _gal *gal, *rand;
  long Ngal, Nrand;
  gal = read_file(ifile,&Ngal);
  rand = read_file(randifile,&Nrand);

  double xmin=1e33;
  double ymin=1e33;
  double zmin=1e33;
  for(long i=0;i<Ngal;i++)
  {
    if(gal[i].x < xmin) xmin = gal[i].x;
    if(gal[i].y < ymin) ymin = gal[i].y;
    if(gal[i].z < zmin) zmin = gal[i].z;
  }
  cout << "Mins: " << xmin << ' '  << ymin << ' ' << zmin << endl;
//Do extra minus to make sure the random x,y,z are also all > 0!
  xmin = xmin-200.;
  ymin = ymin-200.;
  zmin = zmin-200.;
  double wr2sum = 0.;
  double wg2sum = 0.;
  double wg2nzsum=0.;
  double wrsum = 0.;
  double wgsum = 0.;
  double wr2nzsum=0.;
  for(long i=0;i<Ngal;i++)
  {
    double ww,ww2;
    ww = gal[i].weight;
    ww2 = ww*ww;
    gal[i].x = gal[i].x - xmin;
    gal[i].y = gal[i].y - ymin;
    gal[i].z = gal[i].z - zmin;
    wgsum += ww;
    wg2nzsum += ww2*gal[i].nz;
    wg2sum += ww2;
  }
  for(long i=0;i<Nrand;i++)
  {
    double ww,ww2;
    ww = rand[i].weight;
    ww2 = ww*ww;
    rand[i].x = rand[i].x - xmin;
    rand[i].y = rand[i].y - ymin;
    rand[i].z = rand[i].z - zmin;
    wr2sum += ww2;
    wr2nzsum += ww2*rand[i].nz;
    wrsum += ww;
  }
  cout << gal[0].x << ' ' << gal[0].y << ' ' << gal[0].z << ' ' 
       << gal[0].weight << endl;
  cout << rand[0].x << ' ' << rand[0].y << ' ' << rand[0].z << ' ' 
       << rand[0].weight << endl;

  double ***rnrho;
  rnrho = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  double ***nrho;
  nrho = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  cic_mass(rand,Nrand,rnrho,Lbox);
  cic_mass(gal,Ngal,nrho,Lbox);

//Calculate P(k)...
  double ***nrho_in;
  nrho_in = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  double alpha = wgsum/wrsum;
  cout << "alpha = " << alpha << endl;

  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        nrho_in[i][j][k] = nrho[i][j][k]-alpha*rnrho[i][j][k];
      }
    }
  }
  kill_3D_array(nrho,Ngrid,Ngrid);
  kill_3D_array(rnrho,Ngrid,Ngrid);
  cout << "Mean rho: " << calc_mean(nrho_in) << endl;

  double *kk, *pp;
  kk = initialize_vector<double>(nbin);
  pp = initialize_vector<double>(nbin);
  //double norm = wg2nzsum; //Will
  double norm = alpha*wr2nzsum; //Francesco
  //double shot = 1/norm * (wg2sum + alpha*alpha*wr2sum); //Will
  double shot = 1/norm * (alpha*(1.+alpha)*wr2sum); //Francesco
  powspec(nrho_in, nbin, Lbox, norm, shot, kk, pp);
  cout << "shot + norm: " << shot << ' ' << norm << endl;
//  cout << sqrt(w2sum*ng*ng) << ' ' <<  1/ng << ' ' << 1/nr << endl;

  ofstream pkout;
  stringstream ss;
  ss << in << "_Ngrid" << Ngrid << ".pk";
  pkout.open(ss.str().c_str(), ios::out);
  pkout.precision(8);
  for(int i=0;i<nbin;i++) pkout << kk[i] << ' ' << pp[i] << endl;
  pkout.close();
  delete[] kk;
  delete[] pp;

  return 0;
}
