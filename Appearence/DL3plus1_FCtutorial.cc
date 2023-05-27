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

#include "TFile.h"
#include "SBNllminimizer.h"
#include "prob.h"

using namespace sbn;

// define helper functions
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
void printbinedges();

// define some global variables - these are hardcoded to match DL 3+1 joint analysis
const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;
// make output file
std::ofstream chifile;
std::ofstream binsfile;
std::ofstream sensfile;
// set grid points and parameter limits to search for grid search
// dm2 = delta m^2_41 (ev^2)
const int dm2_grdpts(1), ue4_grdpts(5), umu4_grdpts(5);
const double dm2_lowbound(1.89), dm2_hibound(1.89);
const double ue4_lowbound(0.01), ue4_hibound(0.5);
const double umu4_lowbound(0.01), umu4_hibound(0.5);

// this function runs a toy FC method for all of the grid points
// throws 1 universe and returns the minimum -2LLH via a
// minimizer+ -2LLH at the point the universe was thrown from
// A full FC would throw 1k universes at each point and uses this to find the R_crit
// output is a txt file with R, -2LLH bf, -2LLH pt, bf delta m, bf Ue, bf Umu
// in our plotting notebooks in the NueAppearance repository

int main(){

  // --------------------------------------------------
  // initiate more parameters
  // hardcoded - change to your version of the xml
  std::string sbnfithome = "/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star";
  std::string xml = sbnfithome+"/xml/TotalThreePlusOne_full.xml";
  // tag for some saved files
  std::string tag = "DL_full";
  // Number of universes! - set to 10 for this tutorial
  const int nFakeExp(10);
  // --------------------------------------------------


  chifile.open("FC_tutorial.txt", std::ios_base::app);
  binsfile.open("bins_tutorial.txt", std::ios_base::app);
  sensfile.open("sensitivity_tutorial.txt", std::ios_base::app);
  printbinedges();

  // initialize the minimizer object
  SBNllminimizer minimizer( xml );

  // load in saved cvSpec (null oscillation model)
  SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/gen1_oscanalysis/aux/"+tag+"_Bkg.SBNspec.root",xml);
  // Load up the necesary bits and store them in vectors on the stack **
  // first set full osc to 0.0 for the cvspec
  cvSpec.Scale("fullosc",0.0);
  // remove mcerror, since we include it in the systematics matrix
  cvSpec.RemoveMCError();
  std::vector<double> cv_v = cvSpec.collapsed_vector;
  // load in the full systematic matrix - created on fermilab grid
  TFile * fsys = new TFile( "/cluster/tufts/wongjiradlabnu/gen1_oscanalysis/aux/katieversion_bigbins_tot.SBNcovar.root","read");
  TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");

  // create an sbnchi object - used to get cov matrix at given osc + throw universes
  SBNchi cvChi(cvSpec, *covFracSys);
  // save cv covar
  TMatrixD * cvFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
  // fill the object using SBN chiobject
  // uncollapsed martrix has nbins x nbins x number of analysis channels.
  // it needs to be collapsed down into nbins x nbins
  int a = cvChi.FillCollapsedFractionalMatrix(cvFracSys_collapsed);

  // now start the loop through all the grid points....
  // for consistency between scripts, always keep the loops in the same order
  // this will dicate what order grid points are saved to the text files in
  for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
		for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
			for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){
        // now translate to the parameters we are testing
        //mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
        float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
        float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
        float mnu_base = pow(10.,((mi_base+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
        // calculate scaling factors (sin^2(2theta) value) to build the spectra
        float e_app = 4*pow(ue_base,2)*pow(um_base,2);
        float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
        float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));

        // next get the current test model
        NeutrinoModel start_model( mnu_base, uei_base, umui_base );
        minimizer._gen.regenerate_osc( start_model );
        // get the oscillated spectra as the expectation for this grid point
        SBNspec oscSpec =  minimizer.getOscSpectra( mnu_base, ue_base, um_base );
        // again remove SBNfits built in stat error
        oscSpec.RemoveMCError();

        // I'll be generating a new universe around the oscillated spectrum
        // make the chi object and fill cov matrix for this point - the matrix changes with different spectrum contributions
        // inputs SBNspec in, TMatrixT<double> matrix_systematicsin, std::string inxml, bool is_verbose, double random_seed
        // setting the random seed so we always get the same results for the tutorial
        SBNchi TrueChi(oscSpec, *covFracSys, xml, false, 0);
        TrueChi.ReloadCoreSpectrum(&oscSpec);
        fsys->cd();
        TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
        int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);

        // get the motor running with initial cholosky decomposition
        TrueChi.pseudo_from_collapsed = true;
        TrueChi.GeneratePseudoExperiment();

        // Across n fake experiments for this grid point:
        for(int expi = 0; expi < nFakeExp; expi++){
          std::vector<float> fakeData = TrueChi.GeneratePseudoExperiment();

          // get the -2llh at the throw point (pt)
	  TMatrixD cov_pt = GetTotalCov(fakeData, oscSpec, *oscFracSys_collapsed);
          float ptfit = GetLLHFromVector(fakeData, oscSpec, cov_pt, false);

          // now run the minimizer
          std::vector<double> bestfit = minimizer.doFit( fakeData, mnu_base, ue_base, um_base );

          // and save the outputs
          // output is a txt file with R, -2LLH bf, -2LLH pt, bf delta m^2, bf Ue, bf Umu
          chifile<<(ptfit-bestfit[0])<<" "<<bestfit[0]<<" "<<ptfit<<" ";
          chifile<<bestfit[1]<<" "<<bestfit[2]<<" "<<bestfit[3]<<" "
		 <<ue_base<<" "<<um_base<<" "<<mnu_base
		 <<std::endl;

        } //end of fake experiment loop
      } //end of loop over Umu
    } //end of loop over Ue
  } // end of loop over delta m grid


  // next we want to get the info for the sensitivity.
  // (i.e if the data is equal to the null oscillation spectra)
  // We want to calcuate R at each paramter grid point between the CV and the oscillated spectra.
  std::vector<float> data;
  //recast cv doubles to float
  for (int b =0; b<nBins; b++){
    data.push_back(cv_v[b]);
  }

 // now repeat a similar process. loop over all the grid points and find -2LLH at each point
 for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
  for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
    for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){
       // now translate to the parameters we are testing
       //mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
       float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
       float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
       float mnu_base = pow(10.,((mi_base+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
       // calculate scaling factors (sin^2(2theta) value) to build the spectra
       float e_app = 4*pow(ue_base,2)*pow(um_base,2);
       float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
       float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));

       // next get the current test model
       NeutrinoModel start_model( mnu_base, uei_base, umui_base );
       minimizer._gen.regenerate_osc( start_model );
       // get the oscillated spectra as the expectation for this grid point
       SBNspec oscSpec =  minimizer.getOscSpectra( mnu_base, ue_base, um_base );
       // again remove SBNfits built in stat error
       oscSpec.RemoveMCError();
       // make the chi object and fill cov matrix for this point - the matrix changes with different spectrum contributions
       SBNchi TrueChi(oscSpec, *covFracSys);
       TrueChi.ReloadCoreSpectrum(&oscSpec);
       fsys->cd();
       TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
       int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);

       // get the -2llh at the given point
       TMatrixD cov_pt = GetTotalCov(data, oscSpec, *oscFracSys_collapsed);
       float llh = GetLLHFromVector(data, oscSpec, cov_pt, false);

        // and save the outputs
        // output is a txt file with R, -2LLH bf, -2LLH pt, bf delta m^2, bf Ue, bf Umu
       sensfile<<llh<<std::endl;

     } //end of loop over Umu
   } //end of loop over Ue
 } // end of loop over delta m grid
 // now run the minimizer to find the -2llh bf and save to the last time
 std::vector<double> bestfit_sens = minimizer.doFit( data, 0.01, 0.01, 0.01 );
 sensfile<<bestfit_sens[0]<<" "<<bestfit_sens[1]<<" "<<bestfit_sens[2]<<" "<<bestfit_sens[3]<<std::endl;



  return 0;
} // end of main function


// --------------------------------HELPER FUNCTIONS----------------------------

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


void printbinedges(){
	// funtion that prints the bins to the output textfile
	for(int mi = 0; mi <= dm2_grdpts; mi++){
		double mnu = pow(10.,((mi)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		binsfile << pow(mnu,1) << " ";
	}
	binsfile << std::endl;

	for(int uei = 0; uei <= ue4_grdpts; uei++){
		double ue = pow(10.,(uei/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
		binsfile << ue << " ";
	}
	binsfile << std::endl;

	for(int umui = 0; umui <= umu4_grdpts; umui++){
		double umu = pow(10.,(umui/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
		binsfile << umu << " ";
	}
	binsfile << std::endl;
	return;
}//end of print bins function
