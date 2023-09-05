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
#include <algorithm>

// #include <gperftools/profiler.h>

#include "TFile.h"
#include "TRandom3.h"
#include "SBNllminimizer.h"
#include "prob.h"

#include "GridDefinition.h"
#include "KatieOptProcedure.h"

using namespace sbn;
using namespace std::chrono;

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

  // set the input parameter
  int reduced_dimension = atoi(argv[1]);  
  int specific_entry = atoi(argv[2]);
  int start_at_iexp  = -1;
  if ( argc>=4 )
    start_at_iexp = atoi(argv[3]);
  
  std::string reduce_method = "profile";

  // configurations
  std::string sbnfithome = "/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star";
  std::string xml = sbnfithome+"/xml/TotalThreePlusOne_full.xml";  
  
  // Number of universes!
  const int nFakeExp(1000);

  // bins of spectra
  const int nBins_e(12),nBins_mu(19);
  const int nBins = nBins_e+nBins_mu;

  // set grid points: uboone best fit
  const double dm2_data_bestfit = 3.69939;
  const double ue4_data_bestfit = 0.5;
  const double um4_data_bestfit = 0.01;



  GridDefinition griddef;  
  griddef.define_maya_grid();
  int mi_stride   = 20; //  5 pts
  int uei_stride  = 10; //  5 pts
  int umui_stride = 10; //  5 pts
  
  // tag for some save files
  std::string tag = "DL_full";

  // -----------------------------------
  // for testing
  // sample from prior
  // double par_grid1[3] = {0,0,0};
  // if ( reduced_dimension==0 ) {
  //   // this should be log-normal? yes? if our prior is on scale
  //   std::cout << "sample unform between " << log(griddef.dm2_lowbound) << " and " << log(griddef.dm2_hibound) << std::endl;
  //   par_grid1[0] = exp( _rand.Uniform( log(griddef.dm2_lowbound), log(griddef.dm2_hibound)) );
  // }
  // else if (reduced_dimension==1 ) {
  //   par_grid1[1] = _rand.Uniform( griddef.ue4_lowbound, griddef.ue4_hibound);
  // }
  // else if ( reduced_dimension==2 ) {
  //   par_grid1[2] = _rand.Uniform( griddef.umu4_lowbound, griddef.umu4_hibound);
  // }
  
  // std::cout << "Sampled removed parameter: " << par_grid1[ reduced_dimension ] << std::endl;

  // --------------------------------------------------------------------------  

  // now get the parameters we are testing
  // in 2D
  int npars2d = 0;
  int binindex = -1;
  if ( reduced_dimension==0 ) {
    // remove dm2
    npars2d = griddef.ue4_grdpts*griddef.umu4_grdpts;
    if ( specific_entry>=npars2d ) {
      std::cout << "specific entry is out of bounds: " << specific_entry << " vs. npars2d=" << npars2d << std::endl;
      return 0;
    }
    
    // whats the bin index for this specific entry?
    int kindex = 0;
    for(int uei_base = 0; uei_base < griddef.ue4_grdpts; uei_base += 1 ){
      for(int umui_base = 0; umui_base < griddef.umu4_grdpts; umui_base += 1 ) {
	int b = griddef.binindex_to_kindex(0,uei_base,umui_base);
	std::cout << kindex << "==" << specific_entry << " " << b << " " << uei_base << " " << umui_base << std::endl;	  
	
	if ( kindex==specific_entry ) {
	  binindex = griddef.binindex_to_kindex(0,uei_base,umui_base);
	  std::cout << "MATCH: " << kindex << " " << binindex << " " << uei_base << " " << umui_base << std::endl;	  
	}
	kindex++;
	if ( binindex>=0 )
	  break;
      }
      if ( binindex>=0 )
	break;
    }
    std::cout << "kindex: " << kindex << std::endl;
  }
  else if ( reduced_dimension==1 ) {
    // remove Ue4
    npars2d = griddef.dm2_grdpts*griddef.umu4_grdpts;
    if ( specific_entry>=npars2d ) {
      std::cout << "specific entry is out of bounds: " << specific_entry << " vs. npars2d=" << npars2d << std::endl;
      return 0;
    }
    
    // whats the bin index for this specific entry?
    int kindex = 0;
    for(int mi_base = 0; mi_base < griddef.dm2_grdpts; mi_base += 1 ){
      for(int umui_base = 0; umui_base < griddef.umu4_grdpts; umui_base += 1 ) {
	int b = griddef.binindex_to_kindex(mi_base,0,umui_base);
	
	if ( kindex==specific_entry ) {
	  binindex = griddef.binindex_to_kindex(mi_base,0,umui_base);	  
	  break;
	}
	if ( binindex>=0 )
	  break;		
	kindex++;	
      }
      if ( binindex>=0 )
	break;
    }
    std::cout << "kindex: " << kindex << " --> " << binindex << std::endl;
  }
  else if ( reduced_dimension==2 ) {
    // remove Um4
    npars2d = griddef.dm2_grdpts*griddef.ue4_grdpts;
    // whats the bin index for this specific entry?
    int kindex = 0;
    for(int mi_base = 0; mi_base < griddef.dm2_grdpts; mi_base += 1 ){
      for(int uei_base = 0; uei_base < griddef.ue4_grdpts; uei_base += 1 ) {
	if ( kindex==specific_entry ) {
	  binindex = griddef.binindex_to_kindex(mi_base,uei_base,0);
	  break;
	}
	kindex++;
      }
    }        
  }
  
  if ( specific_entry<0 || specific_entry>=npars2d )
    std::runtime_error("Invalide specific entry");

  std::cout << "specific_entry=" << specific_entry << " --> binindex=" << binindex << std::endl;

  if ( binindex<0 ) {
    std::cout << "something wrong in finding the grid point" << std::endl;
    return 0;
  }
  int mi_base   = griddef.binindex_v[binindex][0];
  int uei_base  = griddef.binindex_v[binindex][1];
  int umui_base = griddef.binindex_v[binindex][2];

  float mnu_base = sqrt( griddef.oscpar_v[binindex][0] );
  float ue_base  = griddef.oscpar_v[binindex][1];
  float um_base  = griddef.oscpar_v[binindex][2];

  
  // --------------------------------------------------
  
  // open output files
  std::ofstream chifile;
  std::string textid = ZeroPadNumber(specific_entry, 5);
  chifile.open("chis_seek_"+textid+".txt", std::ios_base::app);

  if ( start_at_iexp<=0 ) {
    chifile << "EXPID" << " "
	    << "DM2" << " " << "Ue4" << " " << "Umu4" << " "
	    << "NLL" << " " << "NLL_chi2" << " " << "NLL_detM" << " "
	    << "NLL-FIT" << " " << "DM2-fit" << " " << "Ue4-fit" << " " << "Um4-fit" << " "
	    << "NLL_fit_chi2" << " " << "NLL_fit_detM"
	    << std::endl;
  }
  
  std::ofstream pseudofile;
  pseudofile.open( "pseudoexp_"+textid+".txt", std::ios_base::app );

  if (start_at_iexp<=0 ) {
    pseudofile << "PSEUDOFITS "
	       << mi_base << " " << uei_base << " " << umui_base << " "
	       << mnu_base << " " << ue_base << " " << um_base
	       << std::endl;
  }
  std::cout << "PSEUDOFITS "
	    << mi_base << " " << uei_base << " " << umui_base << " "
	    << mnu_base << " " << ue_base << " " << um_base
	    << std::endl;

  if (false)
    return 0;  
  
  
  // --------------------------------------------------
  // initiate other parameters
  // hardcoded - change to your version of the xml  
  
  oscfit::KatieOptProcedure analyzer;

  // std::ofstream fitresult_file;
  // fitresult_file.open( "fitresult_"+textid+".txt", std::ios_base::app );

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


  // initialize the minimizer
  // will configure SBNgenerate, which loads the MC events to cache
  SBNllminimizer minimizer( xml );
  std::cout<<"Initialized minimizer"<<std::endl;
  minimizer.use_polar_coords = false;  // never use this    

  // Use the minimizer to regenerate spectrum
  minimizer._gen.regenerate_osc( start_model );
  std::cout << "Full Vector of pre-scale CV spectrum ===========" << std::endl;
  minimizer._gen.spec_central_value.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;
  std::cout << "Full Vector of pre-scale dm2 spectrum ===========" << std::endl;
  minimizer._gen.spec_osc_sinsq.PrintFullVector(true);
  std::cout << "=================================================" << std::endl;

  // print statements to confirm starting spectra
  std::cout << "CV-SPEC: PrintCollapsedVector =======================" << std::endl;
  cvSpec.PrintFullVector();

  float tot_pseudo = 0.0;
  float tot_gridscan = 0.0;
  float tot_simplex = 0.0;
  float tot_migrad = 0.0;

  // Across several fake experiments for this grid point:
  int starting_expi = 0;
  if ( start_at_iexp>0 )
    starting_expi = start_at_iexp;

  TRandom3 _rand(1235+specific_entry+1000*starting_expi);
  
  for(int expi = starting_expi; expi < nFakeExp; expi++){

    // ===================================================
    // START OF PSEUDOEXPERIMENT FIT ROUTINE
    // ===================================================

    auto start_pseudoexp = high_resolution_clock::now();
    std::cout<<"EXPERIMENT NUMBER: "<<expi<<std::endl;

    // for each experiment, we must sample from the prior of the removed variable
    std::vector<double> par_grid = { mnu_base*mnu_base, ue_base, um_base };
    std::vector<int> par_indices = { mi_base, uei_base, umui_base };
    
    // sample from prior
    if ( reduced_dimension==0 ) {
      // this should be log-normal? yes? if our prior is on scale
      std::cout << "sample unform between " << log(griddef.dm2_lowbound) << " and " << log(griddef.dm2_hibound) << std::endl;
      par_grid[0] = exp( _rand.Uniform( log(griddef.dm2_lowbound), log(griddef.dm2_hibound)) ); // is in eV^2
    }
    else if (reduced_dimension==1 ) {
      par_grid[1] = _rand.Uniform( griddef.ue4_lowbound, griddef.ue4_hibound);
    }
    else if ( reduced_dimension==2 ) {
      par_grid[2] = _rand.Uniform( griddef.umu4_lowbound, griddef.umu4_hibound);
    }

    std::cout << "Sampled removed parameter: " << par_grid[ reduced_dimension ] << std::endl;
    
    // get the oscillated spectra as the expectation for this grid point
    // we get a copy of the oscillated spectrum
    SBNspec oscSpec =  minimizer.getOscSpectra( sqrt(par_grid[0]), par_grid[1], par_grid[2] );
    oscSpec.RemoveMCError();

    std::cout << "osc-SPEC: PrintFullVector =======================" << std::endl;
    oscSpec.PrintFullVector(true);

    // Sample a bin count vector
    // make the chi object and fill cov matrix for this point
    SBNchi TrueChi(oscSpec, *covFracSys);
    TrueChi.ReloadCoreSpectrum(&oscSpec);
    TrueChi.InitRandomNumberSeeds( (double)(specific_entry+1000*expi) );
    fsys->cd();
    TMatrixD * oscFracSys_collapsed = (TMatrixD*)fsys->Get("frac_covariance");
    int b = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
    
    // get the motor running with initial cholosky decomposition
    TrueChi.pseudo_from_collapsed = true;
    TrueChi.GeneratePseudoExperiment();

    // sample x, the bin  count vector
    std::vector<float> fakeData = TrueChi.GeneratePseudoExperiment();

    //fitresult_file << "----------------------------" << std::endl;    

    //fitresult_file << std::endl;
    
    // // option to hardcode fakedata spectrum for testing
    // if(false){
    //   fakeData=SetFakeData(fakeData);
    // }

    // get the -2llh at the throw point (pt)
    TMatrixD cov_pt = GetTotalCov(fakeData, oscSpec, *oscFracSys_collapsed);
    
    std::vector<double> oscpar_v = { mnu_base, ue_base, um_base };
    double par_for_nll[3] = { log( par_grid[0] ), par_grid[1], par_grid[2] }; /// { log(dm2), Ue4, Um4 }
    
    // double ptfit = minimizer.negative_likelihood_ratio( par_for_nll );
    // std::vector<double> ptNLL_v = minimizer.llh_components_v;

    std::vector<double> ptNLL_v(3,0);

    if ( reduce_method=="profile" ) {
      analyzer.SKIP_FULL_RES_GRIDFIT = true; 
      minimizer.onlyFitOnePar( reduced_dimension, par_for_nll );
      std::vector<double> pt_result  = analyzer.analyze( fakeData, par_indices,
							 minimizer, griddef,
							 reduced_dimension, "profile" );
      ptNLL_v[0] = pt_result[0];
      ptNLL_v[1] = pt_result[4];
      ptNLL_v[2] = pt_result[5];
      
    }
    else if ( reduce_method=="marginalize" ) {
      ptNLL_v = analyzer.getReducedLikelihood( fakeData, oscpar_v,
					       griddef, minimizer,
					       reduced_dimension,
					       reduce_method,
					       -1, 1000 );
    }
    
    double ptfit = ptNLL_v[0];
    std::cout << "PT_Fit: " << ptfit << std::endl;
    
    //fitresult_file << "ptfit: " << ptfit << std::endl;

    // ===================================================
    // the fit procedures
    // ===================================================

    // run marginalization scheme
    analyzer.SKIP_FULL_RES_GRIDFIT = true;     
    minimizer.fitAllPars();
    std::vector<int> dummy_indices = { 0, 0, 0 };
    std::vector<double> results = analyzer.analyze( fakeData, dummy_indices, minimizer, griddef, -1, "marginalize", 1000 );

    std::cout << "RESULTS: ";
    for ( auto& r : results )
      std::cout << r << " ";
    std::cout << std::endl;
    
    auto stop_pseudoexp = high_resolution_clock::now();
    auto duration_pseudoexp = duration_cast<milliseconds>(stop_pseudoexp - start_pseudoexp);

    // fitresult_file << "time: "
    // 		   << duration_pseudoexp.count() << " "
    // 		   << duration_gridscan.count() << " "
    // 		   << duration_simplex.count() << " "            
    // 		   << duration_migrad.count()
    // 		   << std::endl;

    //tot_pseudo += duration_pseudoexp.count();

    char zelapsed[20];
    sprintf(zelapsed,"%.1f",float(duration_pseudoexp.count())/float(1000.0));

    // save pseudo-experiment info to text file
    pseudofile << "[" << expi << "] ";
    std::cout << "BIN COUNTS for exp[" << expi << "]: ";
    //fitresult_file << "[" << expi << "] ";
    for ( auto& bincount : fakeData ) {
      pseudofile  << bincount << " ";
      std::cout << bincount << " ";
      //fitresult_file << bincount << " ";
    }
    pseudofile << " : " << par_grid[0] << " " << par_grid[1] << " " << par_grid[2];
    pseudofile << " : " << zelapsed << "s";
    pseudofile << std::endl;
    std::cout << std::endl;

    // log results in chi file
    chifile << expi << " " << par_grid[0] << " " << par_grid[1] << " " << par_grid[2] << " ";
    chifile << ptfit << " " << ptNLL_v[1] << " " << ptNLL_v[2] << " ";
    for ( auto& result : results )
      chifile << result << " ";
    chifile << std::endl;
    
  } //end of fake experiment loop

  // fitresult_file << "---------------------------------" << std::endl;
  // fitresult_file << "total times: "
  // 		 << tot_pseudo << " "
  // 		 << tot_gridscan << " "
  // 		 << tot_simplex << " "            
  // 		 << tot_migrad
  // 		 << std::endl;
  
  
  std::cout << "DONE!!" << std::endl;
  //fitresult_file.close();
  
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
