#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <omp.h>
#include <vector>
#include <algorithm>

using namespace std;

#define N_CPU 4

#define L_BOX 3400. //Box size in Mpc/h
#define N_GRID 300 //# of grid cells in each dimension

double lbo2 = L_BOX/2.;
double observer = 0.;

double GRID_SPACE = L_BOX/N_GRID;
long int TOTGRID = N_GRID*N_GRID*N_GRID;

#define RMIN 2. //lower limit of smallest bin
#define RMAX 202. //upper limit of highest bin
#define nR 50 //This is the number of r-bins

/*#define MUMIN -1.
#define MUMAX 1.
#define nMU 40 //This is the number of mu-bins*/

#define MUMIN 0.
#define MUMAX 1.
#define nMU 2 //This is the number of mu-bins

double RMIN2 = RMIN*RMIN, RMAX2=RMAX*RMAX;
double dR = (RMAX-RMIN)/nR;
double dMU = (MUMAX-MUMIN)/nMU;

int CELL_SPAN = ceil(RMAX/GRID_SPACE);

double min(double *x, long int N)
{
  double xmin = 1e33;
  for(long int i=0;i<N;i++)
  {
    if(x[i] < xmin) xmin = x[i];
  }
  return xmin;
}

double max(double *x, long int N)
{
  double xmax = -1e33;
  for(long int i=0;i<N;i++)
  {
    if(x[i] > xmax) xmax = x[i];
  }
  return xmax;
}

class _particle {
  public:
    long int particleID; 
    double x,y,z,weight;
    int which_grid_cell_x(double minx) { return floor( (x-minx)/GRID_SPACE ); }
    int which_grid_cell_y(double miny) { return floor( (y-miny)/GRID_SPACE ); }
    int which_grid_cell_z(double minz) { return floor( (z-minz)/GRID_SPACE ); }
    int grid_cell_number;
};

class _grid {
  public:
    long int nparticles_in_cell;
    long int first_particle_in_cell;
};

class _multipole {
  public:
    void set_vals();
    double DD[nR];
    double DS[nR];
    double SS[nR];
    double RR[nR];
    double xi[nR];
};
void _multipole::set_vals()
{
  for(int i=0;i<nR;i++)
  {
    DD[i] = 0.;
    DS[i] = 0.;
    SS[i] = 0.;
    RR[i] = 0.;
    xi[i] = 0.;
  }
}

/*Assign # of things to be processed by each CPU
   n:		total number of things that need to be processed
   istart:	array of indices corresponding to the object that each cpu 
		will start on 
		NOTE THAT THE LAST ELEMENT OF ISTART == N TO SET WHERE THE
		LAST CPU SHOULD STOP */
int allocate_cpu(long int n, long int istart[N_CPU+1])
{
  for(int i=0;i<N_CPU;i++) istart[i] = n/N_CPU*i;
  istart[N_CPU] = n;

  return 0;
}

//Initialize 1D array
template <class T>
T* initialize_vector(int dim)
{

  T *x;
  x = new T[dim];

  for(int i=0;i<dim;i++) x[i] = 0;

  return x;
}

//Initialize a 3D array
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

double treat_periodic_boundaries(double xin)
{
  double x = xin;

//This line is purely to guard against numerics:
//Sometimes p.z + v.z ~ 0 but slightly negative, without this line, it will add
//L_BOX to p.z and return L_BOX which messes up the gridding (all p.z must be <
//L_BOX!).
    if( x < 1e-4 && x > -1e-4 ) x = 0.;

    if(x < 0.) x += L_BOX;
    if(x >= L_BOX) x -= L_BOX;

    return x;
}

/*Read in particle catalogues
   in:		file name of file containing particle catalogue filenames
   all:		particle objects
   tot:		total number of particles */
int read_particle(string in, vector<_particle> *all, long int *tot)
{
  ifstream hfiles;
  vector<string> hinfile;
  string line;
  long int nfiles,ficpu[N_CPU+1];

  ifstream infile;
  string linel;

  in = in+".txt";
  cout << "Reading " << in << endl;
  infile.open( in.c_str() );

  while( getline(infile,linel) && linel[0] != '#')
  {
    _particle ih;
    stringstream ss(linel);
    ss >> ih.x >> ih.y >> ih.z >> ih.weight;
    ih.z += observer;
  /*  ih.x = treat_periodic_boundaries(ih.x);
    ih.y = treat_periodic_boundaries(ih.y);
    ih.z = treat_periodic_boundaries(ih.z);*/

//    if(rand()/(1.*RAND_MAX) < 0.01) all->push_back(ih);
    all->push_back(ih);
  }
  infile.close();

  *tot = all->size();

  return 0;
}

/*Function for sorting by grid cell*/
int sort_by_grid_number(_particle a, _particle b)
{
  if(a.grid_cell_number < b.grid_cell_number) return 1;
  else return 0;
}

int grid_particles(vector<_particle> *all, long int tot, vector<_grid> *gall,
                   double *mins)
{
  double xmin = mins[0];
  double ymin = mins[1];
  double zmin = mins[2];
  cout << xmin << ' ' << ymin << ' ' << zmin << endl;

  #pragma omp parallel for schedule(dynamic)
  for(long int ip=0;ip<tot;ip++)
  {
    long int i = all->at(ip).which_grid_cell_x(xmin);
    long int j = all->at(ip).which_grid_cell_y(ymin);
    long int k = all->at(ip).which_grid_cell_z(zmin);

    all->at(ip).grid_cell_number = i*N_GRID*N_GRID+j*N_GRID+k;
  }

  sort(all->begin(),all->end(),sort_by_grid_number);

  long int ip=0;
  for(long int cell=0;cell<TOTGRID;cell++)
  {
    _grid tg;
    tg.first_particle_in_cell = ip;
    tg.nparticles_in_cell = 0;
    while(ip < tot && all->at(ip).grid_cell_number == cell)
    {
      tg.nparticles_in_cell++;
      ip++;
    }
    if(tg.nparticles_in_cell == 0) tg.first_particle_in_cell = -1;
    gall->push_back(tg);
//cout << cell << ' ' << gall->at(cell).nparticles_in_cell << ' ' << gall->at(cell).first_particle_in_cell << endl;
  }

  return 0;
}

/*Figure out which cells (around a central cell) we actually need to consider.
If the minimum distance between any corner of a grid cell and any corner of 
the central cell is greater than RMAX, then we don't need to consider it!*/
int skip_cells(bool ***skip)
{
  double cell_size = GRID_SPACE;
  for(int i=-CELL_SPAN;i<CELL_SPAN+1;i++)
  {
    for(int j=-CELL_SPAN;j<CELL_SPAN+1;j++)
    {
      for(int k=-CELL_SPAN;k<CELL_SPAN+1;k++)
      {
        double minsep = 1e33;
        for(int ii=0;ii<=1;ii++){
        for(int jj=0;jj<=1;jj++){
        for(int kk=0;kk<=1;kk++)
        {
          double cor_x = (i+ii)*cell_size;
          double cor_y = (j+jj)*cell_size;
          double cor_z = (k+kk)*cell_size;
          for(int iii=0;iii<=1;iii++){
          for(int jjj=0;jjj<=1;jjj++){
          for(int kkk=0;kkk<=1;kkk++)
          {
            double c_x = iii*cell_size;
            double c_y = jjj*cell_size;
            double c_z = kkk*cell_size;
            double dcx = cor_x - c_x;
            double dcy = cor_y - c_y;
            double dcz = cor_z - c_z;
            double dc = sqrt(dcx*dcx + dcy*dcy + dcz*dcz);
            if(dc < minsep) minsep = dc;
          }}}
        }}}
        if(minsep > RMAX) skip[i+CELL_SPAN][j+CELL_SPAN][k+CELL_SPAN] = 1;
      }
    }
  }

  return 0;
}

//flag = 0 => DD, flag = 1 => DS, flag = 2 => SS, flag = 3 => RR
int count_pairs(vector<_particle> all_a, vector<_particle> all_b,
                long int tot_a, long int tot_b,
                vector<_grid> gall_a, vector<_grid> gall_b,
                double bincount[nR][nMU], bool ***skip, int flag,
                _multipole *xi0, _multipole *xi2, double rat)
{
  cout << "Counting pairs" << endl;
  double bincount_cpu[N_CPU][nR][nMU] = {0};
  double L0count_cpu[N_CPU][nR] = {0};
  double L2count_cpu[N_CPU][nR] = {0};
  double rat2 = rat*rat;

  #pragma omp parallel for schedule(dynamic)
//Loop over grid cells
  for(long int cell=0;cell<TOTGRID;cell++)
  {
    long int cadex = cell;
    if(cadex%1000000==0) cout << "." << flush;
    long int kc = cell%N_GRID;
    long int jc = ((cell-kc)/N_GRID)%N_GRID;
    long int ic = (cell-jc*N_GRID-kc)/(N_GRID*N_GRID);

    long int ipan = gall_a.at(cadex).nparticles_in_cell;
//SKIP grid cell if it's empty
    if(ipan == 0) continue;
    long int ipaf = gall_a.at(cadex).first_particle_in_cell;

    int tid = omp_get_thread_num();
    long int imin = ic-CELL_SPAN;
    long int imax = ic+CELL_SPAN+1;
    long int jmin = jc-CELL_SPAN;
    long int jmax = jc+CELL_SPAN+1;
    long int kmin = kc-CELL_SPAN;
    long int kmax = kc+CELL_SPAN+1;

//Loop over particles in grid cell
    for(long int ipa=ipaf;ipa<ipaf+ipan;ipa++)
    {
      double ax = all_a.at(ipa).x;
      double ay = all_a.at(ipa).y;
      double az = all_a.at(ipa).z;

//Loop over other grid cells of interest
      for(int i=imin;i<imax;i++){
      for(int j=jmin;j<jmax;j++){
      for(int k=kmin;k<kmax;k++)
      {
//Skip cells that are outside RMAX
        if(skip[i-imin][j-jmin][k-kmin]==1) continue;

/*In terms of skipping any double counting, there is a lot to consider. In
particular, the differences between DD/RR and DR are essential to get straight
in one's head.
First off, we actually are _SUPPOSED TO DOUBLE COUNT_!! For example, for DR,
we need to count 1) D on R, then 2) R on D. Same for DD and RR, but 1) and 2)
are the same in these cases, so we can skip double counting and just
multiply the counts by 2 at the end.
For DR, we can use the fact that D on R and R on D will give the same r
separation but the mu separation will have opposite sign. Since we only count
positive mu, we can just consider the absolute value of mu -- this way, we
can just do D on R or R on D (no need to do both!).*/
        if(i < 0 || j < 0 || k < 0 || i >= N_GRID || j >= N_GRID || k >= N_GRID) continue;
        long int cbdex = (i*N_GRID+j)*N_GRID+k;

//Don't double count for DD and RR (_ONLY_ when we are not straddling a 
//boundary --- the flag != 1 part takes care of this)
        if(cbdex < cadex && flag != 1) continue;

//SKIP all empty grid cells
        long int ipbn = gall_b.at(cbdex).nparticles_in_cell;
        if(ipbn == 0) continue;

        long int ipbf = gall_b.at(cbdex).first_particle_in_cell;

//If the cell we're looking at now is the central cell, then we can avoid
//double counting particles by only looping over particles with index > 
//current particle (only for DD and RR, not DR!!)
        long int sd = ipbf;
        if(cbdex == cadex && flag != 1) {sd = ipa+1;}

//Loop over particles in grid cells of interest
        for(long int ipb=sd;ipb<ipbf+ipbn;ipb++)
        {
          double bx = all_b.at(ipb).x;
          double by = all_b.at(ipb).y;
          double bz = all_b.at(ipb).z;
          double abw = all_b.at(ipb).weight*all_a.at(ipa).weight;

          double dx = bx - ax;
          double dy = by - ay;
          double dz = bz - az;

/*          if(dx > lbo2) dx = dx - L_BOX;
          else if(dx < -lbo2) dx = dx + L_BOX;
          if(dy > lbo2) dy = dy - L_BOX;
          else if(dy < -lbo2) dy = dy + L_BOX;
          if(dz > lbo2) dz = dz - L_BOX;
          else if(dz < -lbo2) dz = dz + L_BOX;*/
          double r2 = dx*dx + dy*dy + dz*dz;
          if(r2 < RMIN2 || r2 > RMAX2 || r2 == 0) continue;
          double rsep = sqrt(r2);

          double midx = ax + 0.5*dx; //(ax+bx)/2
          double midy = ay + 0.5*dy;
          double midz = az + 0.5*dz;// + observer;

          double mid2 = midx*midx + midy*midy + midz*midz;
          double midsep = sqrt(mid2); 
          double musep = (midx*dx + midy*dy + midz*dz)/midsep/rsep;
          musep = fabs(musep);
          if(musep < MUMIN || musep > MUMAX) continue;

          int rbin_n = floor((rsep-RMIN)/dR);
          int mubin_n = floor((musep-MUMIN)/dMU);
          if(rbin_n == nR) rbin_n -= 1;
          if(mubin_n == nMU) mubin_n -= 1;

/*cout << cadex << ' ' << cbdex << ' ' << endl;
cout << ipa << ' ' << ipb << ' ' << endl;
cout << rsep << ' ' << rbin_n << ' ' << musep << ' ' << mubin_n << endl;*/
          bincount_cpu[tid][rbin_n][mubin_n] += abw;
          L0count_cpu[tid][rbin_n] += abw;
          L2count_cpu[tid][rbin_n] += 2.5*(3*musep*musep-1.)*abw;
        }
      }}}
    }
  }
  cout << endl;

  for(int i=0;i<N_CPU;i++)
  {
    for(int j=0;j<nR;j++)
    {
      if(flag == 0)
      {
        xi0->DD[j] += L0count_cpu[i][j];
        xi2->DD[j] += L2count_cpu[i][j];
      } 
      else if(flag == 1)
      {
        xi0->DS[j] += L0count_cpu[i][j]*rat;
        xi2->DS[j] += L2count_cpu[i][j]*rat;
      }
      else if(flag == 2)
      {
        xi0->SS[j] += L0count_cpu[i][j]*rat2;
        xi2->SS[j] += L2count_cpu[i][j]*rat2;
      }
      else
      {
        xi0->RR[j] += L0count_cpu[i][j]*rat2;
        xi2->RR[j] += L2count_cpu[i][j]*rat2;
      }
      for(int k=0;k<nMU;k++)bincount[j][k] += bincount_cpu[i][j][k];
    }
  }

  return 0;
}

int initialize_rmu(double rb[nR+1], double mub[nMU+1], 
                   double r[nR], double mu[nMU])
{
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<nR+1;i++) rb[i] = dR*i + RMIN;
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<nMU+1;i++) mub[i] = dMU*i + MUMIN;

  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<nR;i++) r[i] = (rb[i]+rb[i+1])/2.;
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<nMU;i++) mu[i] = (mub[i]+mub[i+1])/2.;

  return 0;
}

int calc_xi(_multipole *xi0, _multipole *xi2)
{
  for(int i=0;i<nR;i++)
  {
/*Since for DD and RR, I avoid double counting, I actually need to multiply
both of them by 2 to get the true counts. I don't need to do this for DR
because I do all of the counts. 
L-S => xi = (DD-2DR+RR)/RR for the true DD, DR and RR counts.
But since, to get the true counts, I need to multiply my DD and RR counts 
by 2, but not my DR counts, this becomes (DD-DR+RR)/RR.*/
    xi0->xi[i] = (xi0->DD[i] - xi0->DS[i] + xi0->SS[i]) / xi0->RR[i];
    xi2->xi[i] = (xi2->DD[i] - xi2->DS[i] + xi2->SS[i]) / xi0->RR[i];
  }
  return 0;
}

int main(int argc, char* argv[])
{
  int ii=1;
  string infile, randfile, randrecfile;
  while(ii < argc)
  {
    string arg = argv[ii];
    if(arg == "-d")
    {
      ii++;
      infile = argv[ii];
      ii++;
    }
    if(arg == "-r")
    {
      ii++;
      randfile = argv[ii];
      ii++;
    }
    if(arg == "-s")
    {
      ii++;
      randrecfile = argv[ii];
      ii++;
    }
  }

  vector<_particle> allparticles, allrandoms, allrandomsSS;
  long int totparticles = 0, totrandoms = 0, totrandomsSS = 0;
  vector<_grid> allgrid, allgridr, allgridrSS;

  cout.precision(10);

  cout << "begin" << endl;
  read_particle(infile,&allparticles,&totparticles);
  cout << "Total number of particles: " << totparticles << endl;
  
  read_particle(randfile,&allrandoms,&totrandoms);
  cout << "Total number of randoms: " << totrandoms << endl;

  read_particle(randrecfile,&allrandomsSS,&totrandomsSS);
  cout << "Total number of shifted randoms: " << totrandomsSS << endl;

//Find minimum x, y, z in particles and randoms for gridding particles
  double *mins,*tt;
  long int totprs = totparticles+totrandoms+totrandomsSS;
  long int totpr = totparticles+totrandoms;
  mins = new double[3];
  tt = new double[totprs];
  cout << "Calculating mins..." << endl;

  for(long int i=0;i<totparticles;i++) tt[i] = allparticles.at(i).x;
  for(long int i=totparticles;i<totpr;i++)
    tt[i] = allrandoms.at(i-totparticles).x;
  for(long int i=totpr;i<totprs;i++)
    tt[i] = allrandomsSS.at(i-totpr).x;
  mins[0] = min(tt,totprs);

  for(long int i=0;i<totparticles;i++) tt[i] = allparticles.at(i).y;
  for(long int i=totparticles;i<totpr;i++)
    tt[i] = allrandoms.at(i-totparticles).y;
  for(long int i=totpr;i<totprs;i++)
    tt[i] = allrandomsSS.at(i-totpr).y;
  mins[1] = min(tt,totprs);

  for(long int i=0;i<totparticles;i++) tt[i] = allparticles.at(i).z;
  for(long int i=totparticles;i<totpr;i++)
    tt[i] = allrandoms.at(i-totparticles).z;
  for(long int i=totpr;i<totprs;i++)
    tt[i] = allrandomsSS.at(i-totpr).z;
  mins[2] = min(tt,totprs);
  delete[] tt;

  grid_particles(&allparticles,totparticles,&allgrid,mins);
  grid_particles(&allrandoms,totrandoms,&allgridr,mins);
  grid_particles(&allrandomsSS,totrandomsSS,&allgridrSS,mins);
  delete[] mins;

  double rbin[nR+1], mubin[nMU+1]; //stores upper/lower limits of bins
  double r[nR], mu[nMU];
  double DD[nR][nMU]={0};
  double DS[nR][nMU]={0};
  double SS[nR][nMU]={0};
  double RR[nR][nMU]={0};
  initialize_rmu(rbin,mubin,r,mu);
  cout << "Initialized" << endl;

  bool ***skip;
  skip = initialize_3D_array<bool>(2*CELL_SPAN+1,2*CELL_SPAN+1,2*CELL_SPAN+1);
  cout << "Figuring out which cells to skip..." << endl;
  skip_cells(skip);

  _multipole *xi0, *xi2;
  xi0 = new _multipole;
  xi2 = new _multipole;
  xi0->set_vals();
  xi2->set_vals();

  double totw=0.;
  for(int i=0;i<totparticles;i++) totw += allparticles.at(i).weight;
  double totwr=0.;
  for(int i=0;i<totrandoms;i++) totwr += allrandoms.at(i).weight;
  double totwrSS=0.;
  for(int i=0;i<totrandomsSS;i++) totwrSS += allrandomsSS.at(i).weight;
  double rat = totw/totwr;
  cout << "Ratio of data weights to random weights: " << rat << endl;
  double ratSS = totw/totwrSS;
  cout << "Ratio of data weights to shifted random weights: " << ratSS << endl;

  count_pairs(allparticles,allparticles,totparticles,totparticles,
              allgrid,allgrid,DD,skip,0,xi0,xi2,rat);
  count_pairs(allparticles,allrandomsSS,totparticles,totrandomsSS,
              allgrid,allgridrSS,DS,skip,1,xi0,xi2,ratSS);
  count_pairs(allrandomsSS,allrandomsSS,totrandomsSS,totrandomsSS,
              allgridrSS,allgridrSS,SS,skip,2,xi0,xi2,ratSS);
  count_pairs(allrandoms,allrandoms,totrandoms,totrandoms,
              allgridr,allgridr,RR,skip,3,xi0,xi2,rat);
  kill_3D_array<bool>(skip,2*CELL_SPAN+1,2*CELL_SPAN+1);

  calc_xi(xi0,xi2);

  for(int i=0;i<nR;i++)
  {
    for(int j=0;j<nMU;j++) cout << r[i] << ' ' << mu[j] << ' ' << DD[i][j] << ' ' << DS[i][j] << ' ' << SS[i][j] << ' ' << RR[i][j] << endl;
  }

  ofstream fout;
  stringstream ss;
  ss << infile << ".xi";
  cout << "Writing " << ss.str() << endl;
  fout.open(ss.str().c_str(), ios::out);
  fout.precision(10);
  for(int i=0;i<nR;i++) 
  {
//Multiply DD and RR counts by 2 to get true counts (since we only did D on D
//and R on R once!)
    fout << r[i] << ' ' << xi0->xi[i] << ' ' << xi2->xi[i] << ' ' 
         << 2.*xi0->DD[i] << ' ' << xi0->DS[i] << ' ' << 2.*xi0->SS[i] << endl;
  }
  cout << "end" << endl;
//  for(int i=0;i<nMU+1;i++) cout << mu[i] << endl;

/*  for(int i=0;i<10;i++)
  {
    cout << allparticles.at(i).x << ' ' << allparticles.at(i).which_grid_cell_x() << endl;
    cout << allparticles.at(i).y << ' ' << allparticles.at(i).which_grid_cell_y() << endl;
    cout << allparticles.at(i).z << ' ' << allparticles.at(i).which_grid_cell_z() << endl;
    cout << allparticles.at(i).grid_cell_number << endl;
  }*/
/*  long int test;
  for(int i=TOTGRID-10;i<TOTGRID;i++)
  {
//      test += allgrid.at(i).nparticles_in_cell;
    cout << allgrid.at(i).first_particle_in_cell << ' ';
    cout << allgrid.at(i).nparticles_in_cell << endl;
  }*/
//  cout << test << endl;
  return 0;
}
