/*!
 * @file multipole-fit.cc
 * @brief Fit multipole moments to a sky map using a series expansion in
 *        spherical harmonics.
 * @author Segev BenZvi
 * @date 30 Apr 2014
 * @version $Id: multipole-fit.cc 19744 2014-05-03 17:34:25Z fiorino $
 */

#include <map-maker/BeamFunctions.h>
#include <map-maker/SmoothMap.h>

#include <data-structures/math/SpecialFunctions.h>

#include <hawcnest/CommandLineConfigurator.h>
#include <hawcnest/Logging.h>
#include <hawcnest/HAWCUnits.h>

#include <healpix_map_fitsio.h>
#include <vec3.h>

#include <TFitter.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace std;
using namespace HAWCUnits;
using namespace SpecialFunctions;

HMap skymap;              // relative intensity
HMap skymapVar;           // relative intensity variance
int minPix;               // minimum pixel ID for looping
int maxPix;               // maximum pixel ID for looping
int lmax;                 // maximum multipole moment to fit
double** ylmTab;          // Lookup table to "memoize" Y_{lm} calculation

static void
chi2func(int &npar, double *gin, double &f, double *alm, int flag)
{
  double sum = 0.;
  double df;
  double sigma2;
  pointing p;

  int pIdx = 0;
  int lIdx = 0;

  // Loop over unmasked pixels
  for (int i = minPix; i <= maxPix; ++i) {
    p = skymap.pix2ang(i);
    df = skymap[i];
    lIdx = 0;

    // Add up spherical harmonic components
    for (int l = 0; l <= lmax; ++l) {
      for (int m = -l; m <= l; ++m) {
        df -= alm[lIdx] * ylmTab[pIdx][lIdx];
        ++lIdx;
      }
    }
    sigma2 = skymapVar[i];
    sum += df*df / sigma2;
    ++pIdx;
  }
  f = sum;
}

HMap Degrade(const HMap orig, int newNSide, bool useAvg = false) {
   HMap map;
   map.SetNside(newNSide, RING);
   map.fill(0.);

   int scale = 1.;
   if (useAvg) {
       scale = map.Npix() / orig.Npix();
   }

   for (int i = 0; i < orig.Npix(); ++i){
       pointing p_i = orig.pix2ang(i);
       int newPixel = map.ang2pix(p_i);
       map[newPixel] += orig[i]*scale;
   }

   return map;

}

int main(int argc, char* argv[])
{
  CommandLineConfigurator cl;
  cl.AddPositionalOption<string>("input", "FITS format map file");
  cl.AddOption<string>("name,n", "hawc", "prefix for map name (e.g., HAWC1Yr");
  cl.AddOption<double>("minDec,m", -26., "min dec in deg (>-90.)");
  cl.AddOption<double>("maxDec,M", 64., "max dec in deg (<=90.)");
  cl.AddOption<double>("smoothing,s", 0., "map smoothing [deg]");
  cl.AddOption<int>("lmax,l", 2, "max multipole moment in fit");
  cl.AddFlag("degrade", "degrade maps to N64");
  if (!cl.ParseCommandLine(argc, argv))
    return 1;

  double minDec = (cl.GetArgument<double>("minDec")) * degree;
  double maxDec = (cl.GetArgument<double>("maxDec")) * degree;
  string input  = cl.GetArgument<string>("input");
  string mapPrefix = cl.GetArgument<string>("name");
  double smoothRadius = cl.GetArgument<double>("smoothing") * degree;
  lmax = cl.GetArgument<int>("lmax");

  // Import keywords from header of first file
  fitshandle in;
  in.open(input);

  double startMJD;     in.get_key("STARTMJD", startMJD);
  double stopMJD;      in.get_key("STOPMJD",  stopMJD);
  double totEvents;    in.get_key("NEVENTS",  totEvents);
  double totEventWgts; in.get_key("SUMWGTS",  totEventWgts);
  double totBkgWgts;   in.get_key("NBKG",     totBkgWgts);
  double totIntHours;  in.get_key("TOTDUR",   totIntHours);
  double avgIntHours;  in.get_key("DURATION", avgIntHours);
  string maptype;      in.get_key("MAPTYPE",  maptype);
  double maxIntHours;  in.get_key("MAXDUR",   maxIntHours);
  double minIntHours;  in.get_key("MINDUR",   minIntHours);

  // Import column names, check to see if the names match what we need
  try {
    in.goto_hdu(2);
    if (in.colname(1) != "relative intensity" &&
        in.colname(2) != "relative intensity error")
    {
      log_fatal_nothrow("Expected FITS column 1 to contain "
                        "'relative intensity,' not '" << in.colname(1) << "'");
    }
  }
  catch (const PlanckError& e) {
    cerr << "Could not access column name in " << input << endl;
    return 1;
  }

  in.close();

  // Import data and background counts, then set up relative intensity
  HMap relInt;    read_Healpix_map_from_fits(input, relInt,    1);
  HMap relIntErr; read_Healpix_map_from_fits(input, relIntErr, 2);
  HMap dat;       read_Healpix_map_from_fits(input, dat,       3);
  HMap bkg;       read_Healpix_map_from_fits(input, bkg,       4);

  // Transform declination limits to pixel limits
  const int npix = dat.Npix();
  pointing pLo(halfpi-minDec, 0.);
  pointing pHi(halfpi-maxDec, 0.);
  minPix = max(0, dat.ang2pix(pHi));
  maxPix = min(npix-1, dat.ang2pix(pLo));

  cout << minPix << " " << maxPix << endl;

  // Fill relative intensity maps and the Y_{lm} lookup table
  skymap.SetNside(relInt.Nside(), relInt.Scheme());
  skymapVar.SetNside(relInt.Nside(), relInt.Scheme());

  const int npar = (lmax + 1)*(lmax + 1);
  const int ndat = maxPix - minPix;
  ylmTab = new double*[ndat];
  int pIdx = 0;
  int lIdx = 0;
  pointing p;

  for (int i = 0; i < bkg.Npix(); ++i) {
    if (i >= minPix && i <= maxPix) {
      // Fill skymap and variance map
      skymap[i] = relInt[i];
      if (isnan(skymap[i]) || isinf(skymap[i]))
        skymap[i] = 0.;
      skymapVar[i] = relIntErr[i]*relIntErr[i];
      if (isnan(skymapVar[i]) || isinf(skymapVar[i]))
        skymapVar[i] = 1e99;

      // Fill Y_{lm} lookup table
      ylmTab[pIdx] = new double[npar];
      lIdx = 0;
      p = skymap.pix2ang(i);
      for (int l = 0; l <= lmax; ++l) {
        for (int m = -l; m <= l; ++m) {
          ylmTab[pIdx][lIdx] = Legendre::ReY(l,m, p.theta,p.phi);
          ++lIdx;
        }
      }
      ++pIdx;
    }
    else {
      skymap[i] = Healpix_undef;
      skymapVar[i] = Healpix_undef;
    }
  }

  // Fit multipole function to the relative intensity map
  TFitter* minimizer = new TFitter(npar);
  minimizer->SetFCN(chi2func);

  int i = 0;
  char parTitle[20];
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      sprintf(parTitle, "a_%d,%d", l, m);
      minimizer->SetParameter(i++, parTitle, 0., 1e-6, -1, 1);
    }
  }

  // Minimize output
  double argList[10];
  argList[0] = 0;
  minimizer->ExecuteCommand("SET PRINTOUT", argList, 1);

  // Set minimizer strategy: 1=standard, 2=improved (slow)
  argList[0] = 2;
  minimizer->ExecuteCommand("SET STR", argList, 1);

  // Set up MIGRAD to iterate up to 10000 times, then print statistics
  argList[0] = 10000;
  minimizer->ExecuteCommand("MIGRAD", argList, 1);

  double amin, edm, errdef;
  int nvpar, nparx;
  minimizer->GetStats(amin, edm, errdef, nvpar, nparx);
  cout << "\nchi2/ndf     = " << amin << "/" << (maxPix-minPix-nvpar)
       << " = " << amin/(maxPix-minPix-nvpar) << endl;
  cout << "\nPr(chi2/ndf) = " << TMath::Prob(amin, maxPix-minPix-nvpar) << endl;

  // Save the fit skymap
  HMap fitmap;
  fitmap.SetNside(skymap.Nside(), skymap.Scheme());

  double beforesum = 0.;
  double aftersum  = 0.;
  pIdx = 0;
  for (int i = 0; i < fitmap.Npix(); ++i) {
    if (i >= minPix && i <= maxPix) {
      // Add up spherical harmonic components
      lIdx = 0;
      fitmap[i] = 0.;
      for (int l = 0; l <= lmax; ++l) {
        for (int m = -l; m <= l; ++m) {
          fitmap[i] += minimizer->GetParameter(lIdx) * ylmTab[pIdx][lIdx];
          ++lIdx;
        }
      }
      ++pIdx;

      // Find the number of signal events needed to get the residual map
      double res = skymap[i] - fitmap[i];
      //dat[i] = (bkg[i] * (res + 1.))*(1-4.98166314424e-5);
      beforesum += dat[i];
      dat[i] = bkg[i] * (res + 1.);
      aftersum += dat[i];
    }
    else {
      fitmap[i] = Healpix_undef;
    }
  }

  double scale = beforesum/aftersum;
  for (int i = 0; i < dat.Npix(); ++i) {
    dat[i] = dat[i]*scale;
  } 

  // Clean up lookup table
  for (int i = 0; i < ndat; ++i)
    delete [] ylmTab[i];
  delete [] ylmTab;

  // Smooth the residual data and background counts and recalculate the
  // relative intensity and significance
  if (cl.HasFlag("degrade")) {
      dat    = Degrade(dat, 64);
      bkg    = Degrade(bkg, 64);
      fitmap = Degrade(fitmap, 64, true);
      minPix = max(0, dat.ang2pix(pHi));
      maxPix = min(npix-1, dat.ang2pix(pLo));
  }

  HMap resmap;
  resmap.SetNside(dat.Nside(), dat.Scheme());

  HMap sigmap;
  sigmap.SetNside(dat.Nside(), dat.Scheme());

/*
  HMap datWSq;
  datWSq.SetNside(dat.Nside(), dat.Scheme());

  HMap bkgWSq;
  bkgWSq.SetNside(dat.Nside(), dat.Scheme());
*/
  double binsize  = .5*sqrt(pi/3.)/(dat.Nside());

  if (smoothRadius > 0.) {
    binsize = smoothRadius; 
/*
    const int lmax = 2*dat.Nside();

    TopHatWSqBeam beamwsq(lmax, dat.Nside(), smoothRadius);
    bkgWSq = SmoothMap::smooth(bkg, beamwsq);
    datWSq = SmoothMap::smooth(dat, beamwsq);

    TopHatBeam beam(lmax, smoothRadius);
    bkg = SmoothMap::smooth(bkg, beam);
    dat = SmoothMap::smooth(dat, beam);
*/
  }

  for (int i = 0; i < bkg.Npix(); ++i) {
    if (i >= minPix && i <= maxPix) {
/*
      double on     = dat[i];
      double off    = bkg[i];
      double onwsq  = datWSq[i];
      double offwsq = bkgWSq[i];
*/

      pointing pt_i = dat.pix2ang(i);
      double dec    =  halfpi - pt_i.theta;

      float off       = 0.;
      float offwsq    = 0.;
      float on        = 0.;
      float onwsq     = 0.;
      /*
      float offraw    = 0.;
      float onROI     = 0.;
      float onwsqROI  = 0.;
      float offROI    = 0.;
      float offwsqROI = 0.;
      */
      vector<int> listOfPixels;
      listOfPixels.clear();
 
      if (smoothRadius != 0) {
 #if HEALPIX_VERSION < 300
        dat.query_disc(pt_i, smoothRadius, listOfPixels);
 #else
        rangeset<int> qdisc;
        dat.query_disc(pt_i, smoothRadius, qdisc);
        for (size_t j = 0; j < qdisc.size(); ++j)
          for (int k = qdisc.ivbegin(j); k < qdisc.ivend(j); ++k)
            listOfPixels.push_back(k);
 #endif
      }
      else
        listOfPixels.push_back(i);

      const size_t size = listOfPixels.size();
      for (size_t j = 0; j < size; ++j){
         int index  = listOfPixels[j];
         on         += dat[index];
         onwsq      += dat[index];//datWSq[index];
         off        += bkg[index];
         offwsq     += bkg[index];//bkgWSq[index];
         /*
         offraw     += rawbkg[index];
         if (roi[index]==1){
           onROI      += dat[index];
           onwsqROI   += datWSq[index];
           offROI     += bkg[index];
           offwsqROI  += bkgWSq[index];
         }
        */
      }

      double alpha  = pi*binsize/ (2.*avgIntHours*(pi/12.)*cos(dec));
      double offcorrected    = off    - alpha*(on    - off);
      double offwsqcorrected = offwsq - alpha*(onwsq - offwsq);
      if (offcorrected<0. || offwsqcorrected<0.) {
        offcorrected    = 0.;
        offwsqcorrected = 0.;
      }
      
      double sigma = 0.;
      if (off > 0.0){
        double scale = 1.;
        if (onwsq + offwsqcorrected > 0.)
          scale = (on+offcorrected/alpha)/(onwsq+offwsqcorrected/alpha);
        double xon   =  on*scale;
        double xoff  =  offcorrected*scale/alpha;
        double Z     = 0.;
        if (xon > 0. && xoff> 0.) { // disallows log(zero)
          double logterm1 = ((1.+alpha)/alpha)*(xon/(xon+xoff));
          double logterm2 = (1.+alpha)*(xoff/(xon+xoff));
          double sqarg = xon*log(logterm1) + xoff*log(logterm2);
          if (sqarg < 0.)
            sqarg = 0.;
          Z = sqrt(2.)*sqrt(sqarg);
          if (offcorrected > xon)
            Z = -Z;
        }   
        sigma = Z;
      }
      sigmap[i] = sigma;
      resmap[i] = (on  - offcorrected) / offcorrected;
    }
    else {
      sigmap[i] = Healpix_undef;
      resmap[i] = Healpix_undef;
    }
  }

  // Write out the fit map
  stringstream ss;
  fitshandle outfit;

  arr<std::string> colname(4);
  colname[0] = "multipole fit";
  colname[1] = "residual intensity";
  colname[2] = "residual data counts";
  colname[3] = "residual significance";

  int w = int(smoothRadius/degree);
  int d = int(1e3 * (smoothRadius/degree - w));
  ss << "!" << mapPrefix << "_N" << fitmap.Nside() << "_S"
     << setw(2) << setfill('0') << w << "p"
     << setw(3) << setfill('0') << d
     << "_multipole_fit.fits.gz";

  outfit.create(ss.str());
  outfit.add_comment("Relative intensity dipole + quadrupole fit map");

  outfit.set_key("STARTMJD", startMJD, "MJD of first event");
  outfit.set_key("STOPMJD",  stopMJD, "MJD of last event");
  outfit.set_key("NEVENTS",  totEvents, "Number of events in map");
  outfit.set_key("SUMWGTS",  totEventWgts, "Summed weights of events in map");
  outfit.set_key("NBKG",     totBkgWgts, "Summed weight of background in map");
  outfit.set_key("SMOOTH",   smoothRadius/degree, "Smoothing radius [deg]");
  outfit.set_key("TOTDUR",   totIntHours, "Total integration time [hr]");
  outfit.set_key("DURATION", avgIntHours, "Average integration time [hr]");
  outfit.set_key("MAPTYPE",  maptype, "e.g., Skymap, Moonmap");
  outfit.set_key("MINDUR",   minIntHours, "Min integration time [hr]");
  outfit.set_key("MAXDUR",   maxIntHours, "Max integration time [hr]");

  size_t idx = input.find_last_of("/");
  idx = (idx == string::npos ? 0 : idx+1);
  outfit.set_key("FITFILE", input.substr(idx), "Map file used for fit");
  outfit.set_key("CHI2", amin, "Fit chi-square");
  outfit.set_key("NDOF", maxPix-minPix-nvpar, "Fit degrees of freedom");
  outfit.set_key("CHI2PROB", TMath::Prob(amin, maxPix-minPix-nvpar),
                 "Fit chi-square probability");

  for (int i = 0; i < minimizer->GetNumberFreeParameters(); ++i) {
    string parname(minimizer->GetParName(i));
    transform(parname.begin(), parname.end(), parname.begin(), ::toupper);
    ss.str("");
    ss << "Fit parameter " << i << endl;
    outfit.set_key(parname, minimizer->GetParameter(i), ss.str()); 
  }

  prepare_Healpix_fitsmap(outfit, fitmap, PLANCK_FLOAT32, colname);
  outfit.write_column(1, fitmap.Map());
  outfit.write_column(2, resmap.Map());
  outfit.write_column(3, dat.Map());
  outfit.write_column(4, sigmap.Map());
  outfit.close();

  cout << endl;
  for (int i = 0; i < minimizer->GetNumberFreeParameters(); ++i) {
    cout << setw(3)  << minimizer->GetParName(i)   << " = "
         << setw(12) << minimizer->GetParameter(i) << " +- "
         << setw(12) << minimizer->GetParError(i)
         << endl;
  }

  return 0;
}
