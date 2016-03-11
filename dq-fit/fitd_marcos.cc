#include <healpix-cxx/healpix/healpix_map.h>
#include <healpix-cxx/healpix/healpix_map_fitsio.h>
#include <healpix-cxx/cxxsupport/fitshandle.h>
#include <healpix-cxx/cxxsupport/vec3.h>

#include <TMinuit.h>

#include <cmath>
#include <iostream>
#include <iomanip>

typedef Healpix_Map<double> HMap;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::sqrt;

const double degree = 4*atan(1.) / 180.;

static HMap map;
static HMap mapVar;

static void
chi2(int& npar, double* gin, double& f, double* par, int iflag)
{
  double m  = par[0]; // Monopole
  double p1 = par[1]; // Dipole
  double p2 = par[2];
  double p3 = par[3];

  double sum = 0;
  double df;
  double sigma2;
  vec3 v;

  const double minZ = sin(-90*degree);
  const double maxZ = sin(-25*degree);

  for (int i = 0; i < map.Npix(); ++i) {
    v = map.pix2vec(i);
    // Mask out the northern hemisphere
    if (v.z >= minZ && v.z <= maxZ) {
      df = map[i] - m -(p1*v.x + p2*v.y + p3*v.z);

      sigma2 = mapVar[i];
      sum += df*df / sigma2;
    }
  }
  f = sum;
}

int main(int argc, char* argv[])
{
  if (argc != 3) {
    cerr << "\nUsage: " << argv[0] << " [bkg.FITS] [dat.FITS]\n\n";
    return 1;
  }

  HMap bkg;
  HMap dat;

  // Read in the relative intensity map
  read_Healpix_map_from_fits(argv[1], bkg);
  read_Healpix_map_from_fits(argv[2], dat);

  map.SetNside(bkg.Nside(), bkg.Scheme());
  mapVar.SetNside(bkg.Nside(), bkg.Scheme());

  HMap mapRes;
  mapRes.SetNside(bkg.Nside(), bkg.Scheme());

  HMap dipmap;
  dipmap.SetNside(bkg.Nside(), bkg.Scheme());

  HMap dqmap;
  dqmap.SetNside(bkg.Nside(), bkg.Scheme());

  double Nb;
  double Nd;
  double minZ = sin(-90*degree);
  double maxZ = sin(-25*degree);
  const double alpha = 1./20.;
  vec3 v;

  // Fill relative intensity maps
  for (int i = 0; i < bkg.Npix(); ++i) {
    v = map.pix2vec(i);
    if (v.z >= minZ && v.z <= maxZ) {
      Nb = bkg[i];
      Nd = dat[i];
      map[i] = (Nd - Nb) / Nb;
      mapVar[i] = Nd*(Nb + alpha*Nd) / (Nb*Nb*Nb);
      if (map[i] != map[i])
        map[i] = 0.;
      if (mapVar[i] != mapVar[i])
        mapVar[i] = 1e99;
    }
    else {
      // Marcos' original values
      //map[i] = 0;
      //mapVar[i] = 0;
      // Segev's values
      map[i] = Healpix_undef;
      mapVar[i] = Healpix_undef;
    }
  }

  // Initialize the MINUIT fit
  TMinuit minuit(4);
  minuit.SetFCN(chi2);

  double argList[10];
  argList[0] = 1;
  int iErrFlag = 0;

  // Set parameter start values, mimizer step size, and min/max limits
  double vstart[4] = { 0, 0, 0, 0 };
  double step[4] = { 1e-6, 1e-6, 1e-6, 1e-6 };
  double bmin[4] = { -1e6, -1e6, -1e6, -1e6 };
  double bmax[4] = {  1e6,  1e6,  1e6, 1e6 };
  const char* name[4] = { "m", "p1", "p2", "p3" };
  for (int i = 0; i < 4; ++i)
    minuit.mnparm(i, name[i], vstart[i], step[i], bmin[i], bmax[i], iErrFlag);

  // Fix a parameter for the fit
  // argList[0] = 0; // Change to parameter number;
  // minuit.mnexcm("FIX", argList, 1, iErrFlag);
  // argList[0] = 8;
  // minuit.mnexcm("FIX", argList, 1, iErrFlag);
  // argList[0] = 4;
  // minuit.mnexcm("FIX", argList, 1, iErrFlag);

  // Set minimization strategy (1=standard; 2=improved, but slow)
  argList[0] = 1;
  minuit.mnexcm("SET STR", argList, 1, iErrFlag);

  // Set up MIGRAD to iterate up to 1000 times
  argList[0] = 100000;
  minuit.mnexcm("MIGRAD", argList, 1, iErrFlag);

  // Scan about the global minimum
  // argList[0] = 0;
  // minuit.mnexcm("SCAN", argList, 1, iErrFlag);

  double x, dx;
  const double scale = 1e3;
  cout << "\nFit Results (x" << 1./scale << "):\n";
  cout.precision(3);
  for (int i = 0; i < 4; ++i) {
    minuit.GetParameter(i, x, dx);
    cout << std::fixed
         << setw(2) << name[i] << " = "
         << setw(6) << x * scale << " +/- "
         << setw(5) << dx * scale
         << endl;
  }


  double m, px, py, pz, dpx, dpy, dpz;
  minuit.GetParameter(0, m, dx);
  minuit.GetParameter(1, px, dpx);
  minuit.GetParameter(2, py, dpy);
  minuit.GetParameter(3, pz, dpz);

  double p = sqrt(px*px + py*py + pz*pz);
  double dp = sqrt((px*px*dpx*dpx + py*py*dpy*dpy + pz*pz*dpz*dpz)/(p*p));

  cout << "dipole: " << p * 1e4 << " +/- " << dp * 1e4 << endl;

  minZ = sin(-90*degree);
  maxZ = sin(-25*degree);

  // Save the residual and dipole maps
  for (int i = 0; i < mapRes.Npix(); ++i) {
    v = mapRes.pix2vec(i);
    // Mask out the northern hemisphere
    mapRes[i] = 0.;
    if (v.z >= minZ && v.z <= maxZ) {
      double dipole = px*v.x + py*v.y + pz*v.z;
      mapRes[i] = map[i] - (m + dipole);
      dipmap[i] = m + dipole;

      // These can't do anything, right?
      //map[i];
      //mapRes[i];
    }
    else {
      mapRes[i] = Healpix_undef;
      dipmap[i] = Healpix_undef;
    }
  }

  fitshandle out1;
  out1.create("!relInt.fits");
  write_Healpix_map_to_fits(out1, map, TDOUBLE);

  fitshandle out2;
  out2.create("!residuals.fits");
  write_Healpix_map_to_fits(out2, mapRes, TDOUBLE);

  fitshandle out3;
  out3.create("!dipolefit.fits");
  write_Healpix_map_to_fits(out3, dipmap, TDOUBLE);

  return 0;
}

