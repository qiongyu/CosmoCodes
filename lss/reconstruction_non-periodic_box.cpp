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

int Ngrid = 256;
int Ngrid3 = Ngrid*Ngrid*Ngrid;
double Lbox = 3400;
int nbin = 100;
double sigma = 15.; //smoothing length for reconstruction
double sigma2 = sigma*sigma;
double bias = 2.0;
double rand_per_cell = 500.;
double f = 0.72;
double mins[3];

struct _gal
{
  float x, y, z, weight;
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
      for(int k=0;k<dim3;k++) x[i][j][k] = 0;
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

int ngp_mass(_gal *p, int N, double ***rho)
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

    rho[i1][i2][i3] += 1;
  }
/*  cout << rho[0][0][0] << endl;
  cout << rho[Ncell][Ncell][Ncell] << endl;
  exit(0);*/

  return 0;
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

int greens(double *deltak, double ***deltax, int dim, int cicflag)
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
        double wcic;
        if(cicflag == 1) wcic = sincx*sincy*sincz;
        else wcic = 1.;

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
/*  int cc=0;
  for(int k=0;k<2*(Ngrid/2+1);k+=2){ cout << deltakd[k] << ' ' << deltakd[k+1]<< endl; }*/

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
        if(k==Ngrid-1)count+=2;
      }
    }
  }
  delete[] deltakd;

  return 0;
}

int reconstruction(double ***delta, double ***delta1, double ***delta2, 
                   double ***delta3, int cicflag)
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
        if(k==Ngrid-1)count+=2;
      }
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
/*  int cc=0;
  for(int k=0;k<2*(Ngrid/2+1);k+=2){ cout << delta_in[k] << ' ' << delta_in[k+1]<< endl; }*/

/*  for(int i=0;i<Ngrid;i++)
  {
    int dex = i*Ngrid*Ngrid+i*Ngrid+i;
    cout << i << ' ' << deltax[dex][0] << ' ' << deltax[dex][1] << ' ' 
         << deltay[dex][0] << ' ' << deltay[dex][1] << ' '
         << deltaz[dex][0] << ' ' << deltaz[dex][1] << endl;
  }*/

//Calculate displacements in k-space...
  greens(delta_in,delta1,0,cicflag);
  greens(delta_in,delta2,1,cicflag);
  greens(delta_in,delta3,2,cicflag);
/*  for(int i=0;i<Ngrid;i++){
  cout << i << ' ' << delta1[i][i][i] << ' ' << delta2[i][i][i] << ' ' << delta3[i][i][i] << endl;
  }*/
  delete[] delta_in;

  return 0;
}

int shift_particles(double ***delta1, double ***delta2, double ***delta3, 
                    int Ngal, _gal *p, int zflag,
                    double ***rand_counts)
{
  if(zflag == 1) cout << "Including RSD correction..." << endl;

//Find displacement at each galaxy position, then shift galaxy back
  double lside = Ngrid/Lbox;
//  #pragma omp parallel for schedule(dynamic) 
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
//    if(rand_counts[i1][i2][i3] < 50.) continue;
    int j1 = i1 + 1;
    int j2 = i2 + 1;
    int j3 = i3 + 1;
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

    p[i].x += dx;
    p[i].y += dy;
    p[i].z += dz;

    if(zflag == 1)
    {
//RSD correction:
      //calculate original position of particle
      double tx = p[i].x + mins[0]; 
      double ty = p[i].y + mins[1];
      double tz = p[i].z + mins[2];
      //calculate distance of particle (LOS distance)
      double los2 = tx*tx + ty*ty + tz*tz; 
      //calculate projection of displacements onto LOS (and multiply by f)
      double dlos = f*(dx*tx + dy*ty + dz*tz)/los2;
      //project back onto x, y, z
      double dx_los = dlos*tx;
      double dy_los = dlos*ty;
      double dz_los = dlos*tz;

      p[i].x += dx_los;
      p[i].y += dy_los;
      p[i].z += dz_los;
    }
  }

  return 0;
}

int calc_delta(string in, double ***delta, _gal *p, int N, 
               double ***rand_counts, int cicflag)
{
  double ***rho;
  rho = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid); 
  if(cicflag == 1)
  {
    cout << "Doing CIC..." << endl;
    cic_mass(p, N, rho);
  }
  else
  {
    cout << "Doing NGP..." << endl;
    ngp_mass(p, N, rho);
  }

//Read in counts file from Mariana
  in = in+".txt";
  cout << "Reading " << in << endl;
  ifstream infile;
  infile.open(in.c_str(), ios::in);

//int count1=0, count2=0;
  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        string line;
        getline(infile,line);
        stringstream ss(line);
        double ra, dec;
        ss >> rand_counts[i][j][k] >> ra >> dec;
//if( (fabs(rand_counts[i][j][k]) < 1e-5) && (fabs(rho[i][j][k]) > 1e-5) ) count1++;
//if( (fabs(rand_counts[i][j][k]) > 1e-5) && (fabs(rho[i][j][k]) < 1e-8) ) count2++;
      }
    }
//    cout << i << ' ' << rand_counts[i][i][i] << ' ' << rho[i][i][i] << endl;
  }
  infile.close();
//cout << "Rands == 0 but rho > 0: " << count1 << endl;
//cout << "Rands > 0 but rho == 0: " << count2 << endl;
//  exit(2);

  double rhomean=0, voleff=0;

  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        if(rand_counts[i][j][k] > 0)
        {
          voleff+=rand_counts[i][j][k]/rand_per_cell;
          rhomean+=rho[i][j][k];
        }
      }
    }
  }
  rhomean = rhomean/voleff;
  cout << "Effective volume: " << voleff << endl;
  cout << "Mean galaxy rho: " << rhomean << endl;

  for(int i=0;i<Ngrid;i++)
  {
    for(int j=0;j<Ngrid;j++)
    {
      for(int k=0;k<Ngrid;k++)
      {
        if(rand_counts[i][j][k] > 0) 
        {
          double fill = rand_counts[i][j][k]/rand_per_cell;
          delta[i][j][k] = rho[i][j][k]/rhomean-fill;
        }
      }
    }
//    cout << i << ' ' << delta[i][i][i] << endl;
  }
  cout << "Mean delta after filling: " << calc_mean(delta) << endl;
  kill_3D_array<double>(rho,Ngrid,Ngrid);
  
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
    double temp;

    getline(infile,line);
    stringstream ss(line);
    ss >> g[i].x >> g[i].y >> g[i].z >> g[i].weight;
  }
  cout << g[0].x << ' ' << g[0].y << ' ' << g[0].z << ' ' 
       << g[0].weight << endl;

  return g;
}

int main(int argc, char* argv[])
{
  int ii=1;
  string in, randin, countin;
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
    if(arg == "-c")
    {
      ii++;
      countin = argv[ii];
      ii++;
    }
  }
  int cicflag=0;
  cout << "Reconstruction using L_box=" << Lbox << "Mpc/h and " << Ngrid
       << " grid cells" << endl;
      
  string ifile = in+".txt";
  string randifile = randin+".txt";

  _gal *gal, *rand;
  int Ngal, Nrand;
  gal = read_file(ifile,&Ngal);
  rand = read_file(randifile,&Nrand);

// SET MINS!!
  mins[0]=-1775.251343-1000.;
  mins[1]=-1633.584595-100.;
  mins[2]=-109.915657-1000.;
  cout << "Mins: " << mins[0] << ' '  << mins[1] << ' ' << mins[2] << endl;
  for(int i=0;i<Ngal;i++)
  {
    gal[i].x = gal[i].x - mins[0];
    gal[i].y = gal[i].y - mins[1];
    gal[i].z = gal[i].z - mins[2];
  }
  for(int i=0;i<Nrand;i++)
  {
    rand[i].x = rand[i].x - mins[0];
    rand[i].y = rand[i].y - mins[1];
    rand[i].z = rand[i].z - mins[2];
  }
  cout << gal[0].x << ' ' << gal[0].y << ' ' << gal[0].z << ' ' 
       << gal[0].weight << endl;
  cout << rand[0].x << ' ' << rand[0].y << ' ' << rand[0].z << ' ' 
       << rand[0].weight << endl;

  double ***delta;
  delta = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  double ***rand_counts;
  rand_counts = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  calc_delta(countin, delta, gal, Ngal, rand_counts, cicflag);
  
//Run reconstruction...
  double ***delta1, ***delta2, ***delta3;
  delta1 = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  delta2 = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);
  delta3 = initialize_3D_array<double>(Ngrid,Ngrid,Ngrid);

  reconstruction(delta,delta1,delta2,delta3,cicflag);

//Shift galaxies...
  cout << "Doing trilinear interpolation on galaxies..." << endl;
  shift_particles(delta1,delta2,delta3,Ngal,gal,1,rand_counts);
  cout << gal[0].x << ' ' << gal[0].y << ' ' << gal[0].z << ' '
       << gal[0].weight << endl;
//Shift randoms...
  cout << "Doing trilinear interpolation on randoms..." << endl;
  shift_particles(delta1,delta2,delta3,Nrand,rand,0,rand_counts);
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
  fout.open(ss.str().c_str(), ios::out);
  fout.precision(10);
  for(int i=0;i<Ngal;i++) fout << gal[i].x + mins[0] << ' ' 
                               << gal[i].y + mins[1] << ' '
                               << gal[i].z + mins[2] << ' ' 
                               << gal[i].weight << endl;
  fout.close();
  ss.str("");

//Output randoms...
  ofstream rfout;
  ss << in << "_rec_Lbox" << Lbox << "_Ngrid" << Ngrid << "_rand.txt";
  rfout.open(ss.str().c_str(), ios::out);
  rfout.precision(10);
  for(int i=0;i<Nrand;i++) rfout << rand[i].x + mins[0] << ' ' 
                                 << rand[i].y + mins[1] << ' '
                                 << rand[i].z + mins[2] << ' ' 
                                 << rand[i].weight << endl;
  fout.close();
  delete[] gal;
  delete[] rand;

  return 0;
}
