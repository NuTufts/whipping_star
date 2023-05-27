/**
 * This routine is for generating a text file of pseudo experiments, returned as events in the bins.
 * The format of each text line is:
 * [ID] [N bins] [values in N bins] [Ue4 for generation] [Um4 for generation] [dm41 for generation]
 * 
 *  1. loads the central value spectrum
 *  2. adjusts the osc parameters
 *  3. runs generator and dumps info into text file
 *
 * Version notes:
 * 2022/09/04 - tmw - modifed from DL3plus1_FCwregen.cc -- it's a single point generation without the fitting!
 *
 **/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <ctime>
#include <chrono>
// #include <gperftools/profiler.h>

#include "TFile.h"
#include "SBNllminimizer.h"
#include "prob.h"

using namespace sbn;

// define helper functions
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
std::string ZeroPadNumber(int num, int digits=3);
std::vector<float> SetFakeData(std::vector<float> fakeData);

// define some global variables
const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;

// this function runs the FC method for a single grid point
// throws specified number of universes and returns the minimum -2LLH via a
// minimizer+ -2LLH at the point the universe was thrown from
// output is a txt file with two lines for each universe

int main(int argc, char* argv[]){

  // --------------------------------------------------
  // initiate other parameters
  // hardcoded - change to your version of the xml
  std::string sbnfithome = "/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star";
  std::string xml = sbnfithome+"/xml/TotalThreePlusOne_full.xml";
  // bins of spectra
  const int nBins_e(12),nBins_mu(19);
  const int nBins = nBins_e+nBins_mu;

  // set grid points: uboone best fit
  // const double dm2_data_bestfit = 3.69939;
  // const double ue4_data_bestfit = 0.5;
  // const double um4_data_bestfit = 0.01;

  // global best fit
  const double dm2_data_bestfit = 1.32;
  const double ue4_data_bestfit = 0.116;
  const double um4_data_bestfit = 0.135;
  
  
  //const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
  //const int dm2_grdpts(1), ue4_grdpts(25), umu4_grdpts(25);  // for best-fit profiles: dm2 fixed
  const int dm2_grdpts(25), ue4_grdpts(1), umu4_grdpts(25);  // for best-fit profiles: Ue4 fixed
  //const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(1);  // for best-fit profiles: Um4 fixed  
  
  const double dm2_lowbound(0.01), dm2_hibound(100);
  //const double dm2_lowbound(dm2_data_bestfit), dm2_hibound(dm2_data_bestfit); // best-fit value
  //const double ue4_lowbound(0.01), ue4_hibound(0.5);
  const double ue4_lowbound(ue4_data_bestfit), ue4_hibound(ue4_data_bestfit);   // best-fit value
  const double umu4_lowbound(0.01), umu4_hibound(0.5);
  //const double umu4_lowbound(um4_data_bestfit), umu4_hibound(um4_data_bestfit);
  // tag for some save files
  std::string tag = "DL_full";
  // Number of universes!
  const int nFakeExp(1000);
  // --------------------------------------------------

  // For profiling
  // we need to align dm2 grid with data-bestfit. we do this by:
  // (1) find the closet grid edge
  // (2) calculate the log(dm^2) shift to align that closest grid to the data bestfit value
  int idx_nearest_dm2 = -1;
  float offset_logdm = 0;
  if ( dm2_grdpts>1 ) {
    float min_offset = 1e9;
    float bestfit_logdm = TMath::Log10(sqrt(dm2_data_bestfit));
    for (int i=0; i<=dm2_grdpts; i++) {
      float grid_logdm = float(i)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound));
      float offset = bestfit_logdm-grid_logdm;
      if ( fabs(offset)<fabs(min_offset) ) {
	min_offset = offset;
	idx_nearest_dm2 = i;
      }
    }
    offset_logdm = min_offset;
  }
  std::cout << "nearest dm2 gridpoint: " << idx_nearest_dm2 << std::endl;
  std::cout << "offset logdm: " << offset_logdm << std::endl;  
  
  // --------------------------------------------------
  
  // initialize the minimizer
  SBNllminimizer minimizer( xml );
  std::cout<<"Initialized minimizer"<<std::endl;

  
  // load in saved cvSpec - change to your location
  SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_Bkg.SBNspec.root",xml);
  // Load up the necesary bits and store them in vectors on the stack **
  cvSpec.Scale("fullosc",0.0);
  // remove mcerror, since we include it in the systematics
  cvSpec.RemoveMCError();
  std::vector<double> cv_v = cvSpec.collapsed_vector;
  // load in the full systematic matrix - created on fermilab grid
  TFile * fsys = new TFile( "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
  TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");

  // create sbnchi - used to get cov matrix at given osc + throw universes
  SBNchi cvChi(cvSpec, *covFracSys);
  // save cv covar for plotting purposes
  TMatrixD * cvFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
  int a = cvChi.FillCollapsedFractionalMatrix(cvFracSys_collapsed);

  // set the oscillation parameters we want to generate psuedo-experiments from
  float ue_base = 0.0; // no nue disappear or appearance
  float um_base = 0.0; // no numu disappearance
  float mnu_base = TMath::Log10(sqrt(dm2_lowbound)); // doesn't matter if no mixing
  // convert back to from log scale to dm
  mnu_base = pow( 10., mnu_base );
  
  // calculate scaling factors (sin^2(2theta) value)
  float e_app = 4*pow(ue_base,2)*pow(um_base,2);
  float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
  float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));

  // get the current test model
  std::cout << "NU Base Model: m41^2:" <<pow(mnu_base,2) <<" m41:"<<mnu_base<< " ue:" << ue_base << " um:" << um_base << std::endl;
  std::cout<<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;
  NeutrinoModel start_model( mnu_base, ue_base, um_base );
  start_model.Printall();
  std::cout << "mass tag: " << start_model.mass_tag << std::endl;

  minimizer._gen.regenerate_osc( start_model );
  std::cout << "Full Vector of pre-scale CV spectrum ===========" << std::endl;
  minimizer._gen.spec_central_value.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;
  std::cout << "Full Vector of pre-scale dm2 spectrum ===========" << std::endl;
  minimizer._gen.spec_osc_sinsq.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;

  // get the oscillated spectra as the expectation for this grid point
  SBNspec oscSpec =  minimizer.getOscSpectra( mnu_base, ue_base, um_base );
  oscSpec.RemoveMCError();

  // I'll be generating a new universe around the oscillated spectrum
  // make the chi object and fill cov matrix for this point
  SBNchi TrueChi(oscSpec, *covFracSys);
  TrueChi.ReloadCoreSpectrum(&oscSpec);
  fsys->cd();
  TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
  int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);

  // print statements to confirm starting spectra
  std::cout << "CV-SPEC: PrintCollapsedVector =======================" << std::endl;
  cvSpec.PrintFullVector();
  std::cout << "osc-SPEC: PrintFullVector =======================" << std::endl;
  oscSpec.PrintFullVector(true);

  // get the motor running with initial cholosky decomposition
  TrueChi.pseudo_from_collapsed = true;
  TrueChi.GeneratePseudoExperiment();

  std::ofstream outfile("pseudo_experiments.csv",std::ofstream::out);

  // Across several fake experiments for this grid point:
  for(int expi = 0; expi < nFakeExp; expi++){
    std::cout<<"EXPERIMENT NUMBER: "<<expi<<std::endl;
    std::vector<float> fakeData = TrueChi.GeneratePseudoExperiment();

    // start of line

    // pseudo-experiment ID, number of bins,
    outfile << expi  << "," << fakeData.size() << ",";
    // the bin values
    for ( auto const binval : fakeData )
      outfile << binval << ",";
    // osc generating points
    outfile << "Ue4:" << ue_base << ",Um4:" << um_base << ",dm:" << mnu_base << std::endl;
      
      
  } //end of fake experiment loop

  outfile.close();
  
  return 0;
  
} // end of main function

TMatrixD GetTotalCov(const std::vector<float>& obsSpec, const SBNspec& expSpec, const TMatrixD& Mfracsys){
	// function to take the fractional Msys and return totl Msys+Mstat
	// inputs:
	// obsSpec: "data" spectra
	// predSpec: "MC" spectra
	// Mfracsys: fractional (flux+xsec+detvar) covariance matrix
	TMatrixD fullcov(nBins,nBins);
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// first set to zero
			fullcov[i][j] = 0.0;
			// scale to the prediction
			fullcov[i][j] = (Mfracsys)[i][j]*expSpec.collapsed_vector[i]*expSpec.collapsed_vector[j];
			// add in stat errors start with CNP for "data" errors;
			if(i==j){
				if (expSpec.collapsed_vector[i] >0){
					// 1/observed 2/expectation
					fullcov[i][j] += 3.0 / (1.0/obsSpec[i] + 2.0/expSpec.collapsed_vector[i]);
				}
				else {
					fullcov[i][j] += expSpec.collapsed_vector[i]/2.0;
				}
			}
		}//end of first bin loop
	}//end of second bin loop
	return fullcov;
}//end of GetTotalCov


float GetLLHFromVector(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Msys, bool prints){
	// // function to calculate a chi2 (shape + rate)
	// // inputs:
	// // obsSpec: "data" vector
	// // predSpec: "MC" spectra
	// // Mfracsys: total (flux+xsec+detvar) covariance matrix
  float chisqTest = 0;
  // inv cov
  TMatrixD invcov = Msys;
  invcov.Invert();

  // add the chi2-like part
  chisqTest = 0;
  for(int i = 0; i < nBins; i++){
    for(int j = 0; j < nBins; j++){
    if(i==j && false) std::cout<<i<<" "
              <<obsSpec[i]<<" "
              <<expSpec.collapsed_vector[i]<<" "
              <<Msys[i][j]<<" "
              <<((obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]))
              <<std::endl;
    chisqTest += (obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]);
        }
      }
  // now need ln(det(2Pi*M))
  chisqTest += log(Msys.Determinant());

  return chisqTest;
}//end of GetLLHFromSpectra

std::string ZeroPadNumber(int num, int digits)
{
  std::stringstream ss;

  // the number is converted to string with the help of stringstream
  ss << num;
  std::string ret;
  ss >> ret;

  // Append zero chars
  int str_length = ret.length();
  for (int i = 0; i < digits - str_length; i++)
    ret = "0" + ret;
  return ret;
}

std::vector<float> SetFakeData(std::vector<float> fakeData){

  // these are the data values!
  fakeData[0] = 4;
  fakeData[1] = 1;
  fakeData[2] = 1;
  fakeData[3] = 2;
  fakeData[4] = 5;
  fakeData[5] = 3;
  fakeData[6] = 8;
  fakeData[7] = 0;
  fakeData[8] = 1;
  fakeData[9] = 0;
  fakeData[10] = 4;
  fakeData[11] = 2;
  fakeData[12] = 26;
  fakeData[13] = 192;
  fakeData[14] = 276;
  fakeData[15] = 401;
  fakeData[16] = 389;
  fakeData[17] = 463;
  fakeData[18] = 439;
  fakeData[19] = 482;
  fakeData[20] = 395;
  fakeData[21] = 353;
  fakeData[22] = 303;
  fakeData[23] = 233;
  fakeData[24] = 240;
  fakeData[25] = 177;
  fakeData[26] = 118;
  fakeData[27] = 109;
  fakeData[28] = 109;
  fakeData[29] = 85;
  fakeData[30] = 58;
  return fakeData;
}
