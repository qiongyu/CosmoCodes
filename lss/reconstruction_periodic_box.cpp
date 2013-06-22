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

int Ngrid = 512;
int Ngrid3 = Ngrid*Ngrid*Ngrid;
double Lbox = 1500;
double sigma = 15.; //smoothing length for reconstruction
double sigma2 = sigma*sigma;
double bias = 2.0;
double f = 0.5;

struct _gal
{
  float x, y, z, vx, vy, vz, weight;
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
T* initialize_vector(int dim)
{

  T *x;
  x = new T[dim];

  for(int i=0;i<dim;i++) x[i] = 0;

  return x;
}

template <class V>
V*** initialize_3D_array(int dim1, int dim2, int dim3)
{

  V ***x;

  x = new V**[dim1];

  for(int i=0;i<dim1;i++)
  {
    x[i] = new V*[dim2];
    for(int j=0;j<dim2;j++)
    {
      x[i][j] = new V[dim3];
      for(int k=0;k<dim3;k++) x[i][j][k] = 0.0;
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
int kill_3D_array(U ***x, int dim1, int dim2)
{
  for(int i=0;i<dim1;i++)
  {
    for(int j=0;j<dim2;j++)
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
  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++) sum += x[i][j][k];
    }
  }
  sum = sum/Ngrid3;
  return sum;
}

int max_pos(_gal *g, int N)
{
  double xmax = -1e33;
  double ymax = -1e33;
  double zmax = -1e33;
  
  for(int i=0;i<N;i++)
  {
    if(g[i].x > xmax) xmax = g[i].x;
    if(g[i].y > ymax) ymax = g[i].y;
    if(g[i].z > zmax) zmax = g[i].z;
  }
  cout << "Max positions: " << xmax << ' ' << ymax << ' ' << zmax << endl;
  return 0;
}

int min_pos(_gal *g, int N)
{
  double xmin = 1e33;
  double ymin = 1e33;
  double zmin = 1e33;
  
  for(int i=0;i<N;i++)
  {
    if(g[i].x < xmin) xmin = g[i].x;
    if(g[i].y < ymin) ymin = g[i].y;
    if(g[i].z < zmin) zmin = g[i].z;
  }
  cout << "Min positions: " << xmin << ' ' << ymin << ' ' << zmin << endl;
  return 0;
}

int get_num_lines(string infile)
{
  int nline=0;
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

int cic_mass(_gal *p, int N, double ***rho)
{
  double lside = Ngrid/Lbox;

  for(int i=0;i<N;i++)
  {
//Because the box is periodic, the Ngrid+1th grid point is the 0th grid point!
//That's why Ngrid == Ncell!
    double dd1 = p[i].x*lside;
    double dd2 = p[i].y*lside;
    double dd3 = p[i].z*lside;
    int i1 = floor(dd1);
    int i2 = floor(dd2);
    int i3 = floor(dd3);
    int j1 = i1 + 1;
    int j2 = i2 + 1;
    int j3 = i3 + 1;
    dd1 = dd1 - i1;
    dd2 = dd2 - i2;
    dd3 = dd3 - i3;
    double de1 = 1. - dd1;
    double de2 = 1. - dd2;
    double de3 = 1. - dd3;
//If the index equals or exceeds Ngrid, then set it back to the 0th or 1st
//grid point!
    if(i1 == Ngrid){i1 = 0; j1 = 1;}
    if(i2 == Ngrid){i2 = 0; j2 = 1;}
    if(i3 == Ngrid){i3 = 0; j3 = 1;}
    if(j1 == Ngrid) j1 = 0;
    if(j2 == Ngrid) j2 = 0;
    if(j3 == Ngrid) j3 = 0;
/*if( (i1 == 1 && i2 == 1 && i3 == 1) || (j1 == 1 && j2 == 1 && j3 == 1) )
{
  cout << p[i].x << ' ' << p[i].y << ' ' << p[i].z << endl;
  cout << dd1 << ' ' << dd2 << ' ' << dd3 << endl;
  cout << de1 << ' ' << de2 << ' ' << de3 << endl;
  cout << i1 << ' ' << i2 << ' ' << i3 << endl;
  cout << j1 << ' ' << j2 << ' ' << j3 << endl;
}*/
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

//    rho[i1][i2][i3] += 1;
    rho[i1][i2][i3] += de1*de2*de3;
    rho[j1][i2][i3] += dd1*de2*de3;
    rho[i1][j2][i3] += de1*dd2*de3;
    rho[j1][j2][i3] += dd1*dd2*de3;
    rho[i1][i2][j3] += de1*de2*dd3;
    rho[j1][i2][j3] += dd1*de2*dd3;
    rho[i1][j2][j3] += de1*dd2*dd3;
    rho[j1][j2][j3] += dd1*dd2*dd3;
  }
/*  cout << rho[0][0][0] << endl;
  cout << rho[Ncell][Ncell][Ncell] << endl;
  exit(0);*/

  return 0;
}

int greens(double *deltak, double ***deltax, int dim)
{
//Set up grid values in k...(don't need this for fft but need for binning)
  double *kgrid;
  double dk = 2.*M_PI/Lbox;
  kgrid = initialize_fft_grid_vectors(dk);

//Calculate displacements, first in k-space, then ifft to configuration space
  double *deltakd;
  deltakd = initialize_vector<double>(Ngrid*Ngrid*(Ngrid+2));

  deltakd[0] = 0.0;
  deltakd[1] = 0.0;

  int r2csizek=Ngrid/2+1;
  double lcico2;
  lcico2 = 0.5*Lbox/Ngrid;
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<Ngrid;i++)
  {
    double sincx = sinc(kgrid[i]*lcico2);
    sincx *= sincx;
    for(int j=0;j<Ngrid;j++)
    {
      double sincy = sinc(kgrid[j]*lcico2);
      sincy *= sincy;
      for(int k=0;k<r2csizek;k++)
      {
        double sincz = sinc(kgrid[k]*lcico2);
        sincz *= sincz;
//For deconvolving the CIC window function:
        double wcic = sincx*sincy*sincz; 

        int dex = (i*Ngrid + j)*r2csizek + k;
        int redex = 2*dex;
        int imdex = redex+1;
        double re = deltak[redex]/Ngrid3 /wcic;
        double im = deltak[imdex]/Ngrid3 /wcic;

        if(i==0 && j==0 && k==0) continue;

//Calculates displacements in k-space for each dimension...
//Note the smoothing by a Gaussian filter
        double k2 = kgrid[i]*kgrid[i] + kgrid[j]*kgrid[j] + kgrid[k]*kgrid[k];
        double smoothok2 = -exp_x(-sigma2*k2/2.)/k2;

        double kuse;
        if(dim==0) kuse=kgrid[i];
        else if(dim==1) kuse=kgrid[j];
        else kuse=kgrid[k];
        deltakd[redex] = -kuse *im*smoothok2;
        deltakd[imdex] =  kuse *re*smoothok2;
/*cout << i << ' ' << j << ' ' << k << ' ' << endl;
cout << redex << ' ' << imdex << endl;
cout << deltak[redex] << ' ' << deltak[imdex] << endl;
cout << Ngrid3 << endl;
cout << k2 << ' ' << smoothok2 << endl;
cout << kuse << endl;
cout << deltakd[redex] << ' ' << deltakd[imdex] << endl;
exit(2);*/
      }
    }
  }
  delete[] kgrid;
/*for(int i=0;i<Ngrid;i++){ int dex = Ngrid*Ngrid*i + Ngrid*i + i; 
cout << i <<  ' ' << deltakx[dex][0] << ' ' << deltakx[dex][1] << endl;}
exit(2);*/
 
//Inverse FFT to get back the displacements in real space... 
  fftw_init_threads();
  fftw_plan_with_nthreads(4);

  fftw_plan plan_d;
  plan_d = fftw_plan_dft_c2r_3d(Ngrid,Ngrid,Ngrid,(fftw_complex *)deltakd,
                                deltakd,FFTW_ESTIMATE);

  cout << "Doing backwards transform on dim " << dim << "..." << endl;
  fftw_execute(plan_d);
  fftw_destroy_plan(plan_d);
  fftw_cleanup_threads();

  int count=0;

  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        deltax[i][j][k] = deltakd[count];
        count++;
      }
      count+=2;
    }
//    cout << i << ' ' << delta1[i][i][i] << ' '  << delta2[i][i][i] << ' '
//         << delta3[i][i][i] << endl;
  }
  delete[] deltakd;

  return 0;
}

int reconstruction(double ***delta, double ***delta1, double ***delta2, 
                   double ***delta3)
{
//Do fft to get delta_k...
  double *delta_in;
  int r2csize=Ngrid*Ngrid*(Ngrid+2);
  delta_in = initialize_vector<double>(r2csize);
  int count=0;

  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        delta_in[count] = delta[i][j][k]/bias;
        count++;
//For r2c, array padding (2 elements) at the end of each k-row. So skip 2 if
//k = Ngrid-1 (last element of k-row).
      }
      count+=2;
    }
  }
  assert(count == r2csize);

  cout << "Calculating the FFT..." << flush;

  fftw_init_threads();
  fftw_plan_with_nthreads(4);

  fftw_plan plan;
//Do fft in-place!
  plan = fftw_plan_dft_r2c_3d(Ngrid,Ngrid,Ngrid,delta_in,
                               (fftw_complex *)delta_in, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup_threads();

  cout << "Done with FFT...on to calculating displacements" << endl;

/*  for(int i=0;i<Ngrid;i++)
  {
    int dex = i*Ngrid*Ngrid+i*Ngrid+i;
    cout << i << ' ' << deltax[dex][0] << ' ' << deltax[dex][1] << ' ' 
         << deltay[dex][0] << ' ' << deltay[dex][1] << ' '
         << deltaz[dex][0] << ' ' << deltaz[dex][1] << endl;
  }*/

//Calculate displacements in k-space...
  greens(delta_in,delta1,0);
  greens(delta_in,delta2,1);
  greens(delta_in,delta3,2);
  delete[] delta_in;

  return 0;
}

int shift_particles(double ***delta1, double ***delta2, double ***delta3, 
                    int Ngal, _gal *p, int zflag)
{
//Find displacement at each galaxy position, then shift galaxy back
  double lside = Ngrid/Lbox;
  #pragma omp parallel for schedule(dynamic) 
  for(int i=0;i<Ngal;i++)
  {
//Because the box is periodic, the Ngrid+1th grid point is the 0th grid point!
//That's why Ngrid == Ncell!
    double dd1 = p[i].x*lside;
    double dd2 = p[i].y*lside;
    double dd3 = p[i].z*lside;
    int i1 = floor(dd1);
    int i2 = floor(dd2);
    int i3 = floor(dd3);
    int j1 = i1 + 1;
    int j2 = i2 + 1;
    int j3 = i3 + 1;
    dd1 = dd1 - i1;
    dd2 = dd2 - i2;
    dd3 = dd3 - i3;
    double de1 = 1. - dd1;
    double de2 = 1. - dd2;
    double de3 = 1. - dd3;
//If the index equals or exceeds Ngrid, then set it back to the 0th or 1st
//grid point!
    if(i1 == Ngrid){i1 = 0; j1 = 1;}
    if(i2 == Ngrid){i2 = 0; j2 = 1;}
    if(i3 == Ngrid){i3 = 0; j3 = 1;}
    if(j1 == Ngrid) j1 = 0;
    if(j2 == Ngrid) j2 = 0;
    if(j3 == Ngrid) j3 = 0;
    assert(i1 >= 0 && i1 < Ngrid &&
           i2 >= 0 && i2 < Ngrid &&
           i3 >= 0 && i3 < Ngrid &&
           j1 >= 0 && j1 < Ngrid &&
           j2 >= 0 && j2 < Ngrid &&
           j3 >= 0 && j3 < Ngrid);

//do 3D linear interpolation to get dx, dy and dz at position of galaxy
    double d00 = de1*delta1[i1][i2][i3] + dd1*delta1[j1][i2][i3];
    double d10 = de1*delta1[i1][j2][i3] + dd1*delta1[j1][j2][i3];
    double d01 = de1*delta1[i1][i2][j3] + dd1*delta1[j1][i2][j3];
    double d11 = de1*delta1[i1][j2][j3] + dd1*delta1[j1][j2][j3];
    double d0 = d00*de2 + d10*dd2;
    double d1 = d01*de2 + d11*dd2;
    double dx = d0*de3 + d1*dd3;

    d00 = de1*delta2[i1][i2][i3] + dd1*delta2[j1][i2][i3];
    d10 = de1*delta2[i1][j2][i3] + dd1*delta2[j1][j2][i3];
    d01 = de1*delta2[i1][i2][j3] + dd1*delta2[j1][i2][j3];
    d11 = de1*delta2[i1][j2][j3] + dd1*delta2[j1][j2][j3];
    d0 = d00*de2 + d10*dd2;
    d1 = d01*de2 + d11*dd2;
    double dy = d0*de3 + d1*dd3;

    d00 = de1*delta3[i1][i2][i3] + dd1*delta3[j1][i2][i3];
    d10 = de1*delta3[i1][j2][i3] + dd1*delta3[j1][j2][i3];
    d01 = de1*delta3[i1][i2][j3] + dd1*delta3[j1][i2][j3];
    d11 = de1*delta3[i1][j2][j3] + dd1*delta3[j1][j2][j3];
    d0 = d00*de2 + d10*dd2;
    d1 = d01*de2 + d11*dd2;
    double dz = d0*de3 + d1*dd3;

    p[i].x = p[i].x + dx;
    p[i].y = p[i].y + dy;
    p[i].z = p[i].z + dz;
    if(zflag == 1) p[i].z = p[i].z + dz*f;
  }

  return 0;
}

int calc_delta(double ***delta, _gal *p, int N)
{
  double ***rho;
  rho = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid); 
  cic_mass(p, N, rho);
  for(int i=0;i<Ngrid;i++)  cout << i << ' '  << rho[i][i][i] << endl;

  double rhomean;
  rhomean = calc_mean(rho);
  cout << "Mean galaxy rho: " << rhomean << endl;

  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        delta[i][j][k] = rho[i][j][k]/rhomean-1.;
      }
    }
 //   cout << i << ' ' << delta[i][i][i] << ' ' << rho[i][i][i] << endl;
  }
//  cout << rho[0][0][1] << ' ' << rho[0][1][1] << ' ' << rho[1][0][1] << ' ' << rho[1][1][0] << ' ' <<  rho[0][1][0] << ' ' << rho[1][0][0] << endl;
  kill_3D_array<double>(rho,Ngrid,Ngrid);
  
  return 0;
}

double treat_periodic_boundaries(double xin)
{
  double x = xin;

//This line is purely to guard against numerics:
//Sometimes p.z + v.z ~ 0 but slightly negative, without this line, it will add
//Lbox to p.z and return Lbox which messes up the gridding (all p.z must be <
//Lbox!).
  if( x < 1e-4 && x > -1e-4 ) x = 0.;

  if(x < 0.) x += Lbox;
  if(x >= Lbox) x -= Lbox;

  return x;
}


int z_distort(_gal *g, int N)
{
  cout << "Doing RSD..." << endl;

  for(int i=0;i<N;i++){
    g[i].z = g[i].z + g[i].vz;
    g[i].z = treat_periodic_boundaries(g[i].z);
  }

  return 0;
}

_gal* read_file(string fname, int *N)
{
  //Read in galaxies
  ifstream infile;

  cout << "Reading " << fname << endl;
  *N = get_num_lines(fname);
  cout << "Number of objects = " << *N << endl;

  _gal *g;
  g = new _gal[*N];
  infile.open( fname.c_str(), ios::in );

  for(int i=0;i<*N;i++)
  {
    string line;

    getline(infile,line);
    stringstream ss(line);
    ss >> g[i].x >> g[i].y >> g[i].z;
    g[i].weight = 1.;
  }

  return g;
}

int main(int argc, char* argv[])
{
  int ii=1;
  int zflag = 0;
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
    if(arg == "-z")
    {
      ii++;
      zflag = 1;
    }
  }
      
  string ifile = in+".txt";
  string randifile = randin+".txt";

  _gal *gal, *rand;
  int Ngal, Nrand;
  gal = read_file(ifile,&Ngal);
  cout << "Galaxies:" << endl;
  cout << gal[0].x << ' ' << gal[0].y << ' ' << gal[0].z << endl; 
  rand = read_file(randifile,&Nrand);
  cout << "Randoms:" << endl;
  cout << rand[0].x << ' ' << rand[0].y << ' ' << rand[0].z << endl; 

  if(zflag == 1)
  {
    z_distort(gal,Ngal);
    cout << "RSD Galaxies:" << endl;
    cout << gal[0].x << ' ' << gal[0].y << ' ' << gal[0].z << endl;
  }
  max_pos(gal,Ngal);
  min_pos(gal,Ngal);

//Calculate delta at grid points (via CIC and filling in empty cells)...
  double ***delta;
  delta = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  calc_delta(delta, gal, Ngal);
  
//Run reconstruction...
  double ***delta1, ***delta2, ***delta3;
  delta1 = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  delta2 = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  delta3 = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);

  reconstruction(delta,delta1,delta2,delta3);

//Shift galaxies...
  cout << "Doing trilinear interpolation on galaxies..." << endl;
  shift_particles(delta1,delta2,delta3,Ngal,gal,zflag);
  cout << "Shifted Galaxies:" << endl;
  cout << gal[0].x << ' ' << gal[0].y << ' ' << gal[0].z << ' '
       << gal[0].weight << endl;
//Shift randoms...
  cout << "Doing trilinear interpolation on randoms..." << endl;
  shift_particles(delta1,delta2,delta3,Nrand,rand,0);
  cout << "Shifted Randoms:" << endl;
  cout << rand[0].x << ' ' << rand[0].y << ' ' << rand[0].z << ' '
       << rand[0].weight << endl;

  kill_3D_array<double>(delta,Ngrid,Ngrid);
  kill_3D_array<double>(delta1,Ngrid,Ngrid);
  kill_3D_array<double>(delta2,Ngrid,Ngrid);
  kill_3D_array<double>(delta3,Ngrid,Ngrid);

//Output galaxies...
  ofstream fout;
  stringstream ss;
  ss << in << "_rec_Lbox" << Lbox << "_Ngrid" << Ngrid << ".txt";
  cout << "Writing " << ss.str() << endl;
  fout.open(ss.str().c_str(), ios::out);
  fout.precision(10);
  for(int i=0;i<Ngal;i++) fout << gal[i].x << ' ' 
                               << gal[i].y << ' '
                               << gal[i].z << endl;
  fout.close();
  ss.str("");
  delete[] gal;

//Output galaxies...
  ofstream rfout;
  ss << in << "_rec_Lbox" << Lbox << "_Ngrid" << Ngrid << "_rand.txt";
  cout << "Writing " << ss.str() << endl;
  rfout.open(ss.str().c_str(), ios::out);
  rfout.precision(10);
  for(int i=0;i<Nrand;i++) rfout << rand[i].x << ' ' 
                                 << rand[i].y << ' '
                                 << rand[i].z << endl;
  rfout.close();
  delete[] rand;

  return 0;
}
