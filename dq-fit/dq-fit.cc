
#include <healpix-cxx/healpix/healpix_map.h>
#include <healpix-cxx/healpix/healpix_map_fitsio.h>
#include <healpix-cxx/cxxsupport/fitshandle.h>
#include <healpix-cxx/cxxsupport/vec3.h>

#include <TFitter.h>
#include <TMath.h>
#include <TMinuit.h>

#include <boost/format.hpp>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

typedef Healpix_Map<double> HMap;

const double deg2rad = 4*atan(1.) / 180.;

static HMap skymap;
static HMap skymapVar;

// Marcos' script goes from -90 to -25
const double minZ = sin(-88*deg2rad);
const double maxZ = sin(-25*deg2rad);

static void
chi2(int& npar, double* gin, double& f, double* par, int iflag)
{
  double m  = par[0]; // Monopole
  double p1 = par[1]; // Dipole
  double p2 = par[2];
  double p3 = par[3];
  double Q1 = par[4]; // Quadrupole
  double Q2 = par[5];
  double Q3 = par[6];
  double Q4 = par[7];
  double Q5 = par[8];

  double sum = 0;
  double df;
  double sigma2;
  vec3 v;

  for (int i = 0; i < skymap.Npix(); ++i) {
    v = skymap.pix2vec(i);
    // Mask out the northern hemisphere
    if (v.z >= minZ && v.z <= maxZ) {
      df = skymap[i] - m -(p1*v.x + p2*v.y + p3*v.z +
                         Q1 * 0.5*(3*v.z*v.z-1.) +
                         Q2 * 2*v.x*v.z +
                         Q3 * 2*v.y*v.z +
                         Q4 * (v.x*v.x - v.y*v.y) +
                         Q5 * 2*v.x*v.y);

      sigma2 = skymapVar[i];
      sum += df*df / sigma2;
    }
  }
  f = sum;
}

void
OutputFitToTex(const TFitter* minimizer, const double scale)
{
  using boost::format;

  const int nPars = minimizer->GetNumberTotalParameters();
  double corr;

  cout << "\n";
  for (int i = 0; i < nPars; ++i) {
    cout << format("$%s$ & $%6.3f\\pm%.3f$")
            % minimizer->GetParName(i)
            % (scale * minimizer->GetParameter(i))
            % (scale * minimizer->GetParError(i));
    for (int j = 0; j < nPars; ++j) {
      if (j <= i) {
        corr = minimizer->GetCovarianceMatrixElement(i, j);
        corr /= sqrt(minimizer->GetCovarianceMatrixElement(i, i));
        corr /= sqrt(minimizer->GetCovarianceMatrixElement(j, j));
        cout << " & " << format("%6.3f") % corr;
      }
      else {
        cout << " &       ";
      }
    }
    cout << "\\\\" << endl;
    // cout << (i+1<nPars ? " \\\\\n" : "\n");
  }
}

void
OutputFit(const TFitter* minimizer, const double scale, const int ndata)
{
  using boost::format;

  double amin;
  double edm;
  double errdef;
  int nvpar;
  int nparx;
  minimizer->GetStats(amin, edm, errdef, nvpar, nparx);

  cout << format("\nchi2 / ndf = %.1f / %d = %.2g\n\n")
          % amin % (ndata-nvpar)
          % TMath::Prob(amin, ndata-nvpar);

  cout << "Fit values (x" << scale << "):\n";
  for (int i = 0; i < nparx; ++i) {
    cout << format("%s = %6.3f\n")
            % minimizer->GetParName(i)
            % (scale * minimizer->GetParameter(i));
  }

  const double px  = minimizer->GetParameter(1);
  const double dpx = minimizer->GetParError(1);
  const double py  = minimizer->GetParameter(2);
  const double dpy = minimizer->GetParError(2);
  const double pz  = minimizer->GetParameter(3);
  const double dpz = minimizer->GetParError(3);

  const double p = sqrt(px*px + py*py + pz*pz);
  const double dp = sqrt((px*px*dpx*dpx + py*py*dpy*dpy + pz*pz*dpz*dpz)/(p*p));

  cout << "\nDipole amplitude (x" << scale << "): "
       << scale*p << " +/- " << scale*dp << endl;
}

int main(int argc, char* argv[])
{
  if (argc != 4) {
    cerr << "\nUsage: " << argv[0] << " [input.FITS] [outFile] [output] \n\n";
    return 1;
  }

  HMap bkg;
  HMap dat;

  // Read in data and background maps
  read_Healpix_map_from_fits(argv[1], dat, 1);
  read_Healpix_map_from_fits(argv[1], bkg, 2);
  string outFile(argv[2]);
  string output(argv[3]);

  skymap.SetNside(bkg.Nside(), bkg.Scheme());
  skymapVar.SetNside(bkg.Nside(), bkg.Scheme());

  HMap skymapRes;
  skymapRes.SetNside(bkg.Nside(), bkg.Scheme());

  double Nb;
  double Nd;
  const double alpha = 1./20.;
  vec3 v;
  int ndata = 0;

  for (int i = 0; i < bkg.Npix(); ++i) {
    v = skymap.pix2vec(i);
    if (v.z >= minZ && v.z <= maxZ) {
      Nb = bkg[i];
      Nd = dat[i];
      skymap[i] = (Nd - Nb) / Nb;
      skymapVar[i] = Nd*(Nb + alpha*Nd) / (Nb*Nb*Nb);
      if (skymap[i] != skymap[i])
        skymap[i] = 0.;
      if (skymapVar[i] != skymapVar[i])
        skymapVar[i] = 1e99;
      ++ndata;
    }
    else {
      skymap[i] = 0;
      skymapVar[i] = 0;
    }
  }

  // Initialize the MINUIT fit
  TFitter* minimizer = new TFitter(9);
  minimizer->SetFCN(chi2);

  // Reduce output (!!!)
  double argList[10];
  argList[0] = -1;
  minimizer->ExecuteCommand("SET PRINTOUT", argList, 1);

  minimizer->SetParameter(0, "m_0", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(1, "p_x", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(2, "p_y", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(3, "p_z", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(4, "Q_1", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(5, "Q_2", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(6, "Q_3", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(7, "Q_4", 0., 1e-6, -1e6, 1e6);
  minimizer->SetParameter(8, "Q_5", 0., 1e-6, -1e6, 1e6);

  // Suppress the dipole terms
  //minimizer->FixParameter(1);
  //minimizer->FixParameter(2);
  //minimizer->FixParameter(3);

  // Suppress the quadrupole terms
  //minimizer->FixParameter(4);
  //minimizer->FixParameter(5);
  //minimizer->FixParameter(6);
  //minimizer->FixParameter(7);
  //minimizer->FixParameter(8);

  // Set the minimization strategy (1=standard, 2=improved (slow))
  argList[0] = 2;
  minimizer->ExecuteCommand("SET STR", argList, 1);

  // Set up MIGRAD to iterate up to 1000 times
  argList[0] = 1000;
  minimizer->ExecuteCommand("MIGRAD", argList, 1);

  if (output != "false") {
    OutputFitToTex(minimizer, 1e4);
    OutputFit(minimizer, 1e4, ndata);
  }

  // Save the fit and residual skymap
  HMap fitmap;
  fitmap.SetNside(skymap.Nside(), skymap.Scheme());

  const double m0 = minimizer->GetParameter(0);
  const double px = minimizer->GetParameter(1);
  const double py = minimizer->GetParameter(2);
  const double pz = minimizer->GetParameter(3);
  double* Q = NULL;
  Q = new double[6];
  Q[1] = minimizer->GetParameter(4);
  Q[2] = minimizer->GetParameter(5);
  Q[3] = minimizer->GetParameter(6);
  Q[4] = minimizer->GetParameter(7);
  Q[5] = minimizer->GetParameter(8);

  for (int i=0; i<fitmap.Npix(); ++i) {
    v = fitmap.pix2vec(i);
    if (v.z >= minZ && v.z <= maxZ) {
      fitmap[i] = m0 + px*v.x + py*v.y + pz*v.z;
      fitmap[i] += Q[1] * 0.5*(3*v.z*v.z-1.) +
                   Q[2] * 2*v.x*v.z +
                   Q[3] * 2*v.y*v.z +
                   Q[4] * (v.x*v.x - v.y*v.y) +
                   Q[5] * 2*v.x*v.y;
    }
  }

  if (outFile != "false") {

    stringstream sstr;
    sstr.str("");
    fitshandle out;

    sstr << outFile;
    out.create(sstr.str().c_str());
    write_Healpix_map_to_fits(out, fitmap, FITSUTIL<double>::DTYPE);
    out.close();
  }

  return 0;
}
