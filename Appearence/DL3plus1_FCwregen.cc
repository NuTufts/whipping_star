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
#include "SBNllminimizer.h"
#include "prob.h"

#include "GridDefinition.h"

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
  int specific_entry = atoi(argv[1]);

  bool DO_SIMPLEX_SEARCH = true;
  const int NTOP_GRID_TO_FIT = 3;

  // Number of universes!
  const int nFakeExp(1000);

  // --------------------------------------------------
  // initiate other parameters
  // hardcoded - change to your version of the xml
  std::string sbnfithome = "/cluster/tufts/wongjiradlabnu/kmason03/taritreetests/whipping_star";
  std::string xml = sbnfithome+"/xml/TotalThreePlusOne_full.xml";
  // bins of spectra
  const int nBins_e(12),nBins_mu(19);
  const int nBins = nBins_e+nBins_mu;

  // set grid points: uboone best fit
  const double dm2_data_bestfit = 3.69939;
  const double ue4_data_bestfit = 0.5;
  const double um4_data_bestfit = 0.01;

  // global best fit
  // const double dm2_data_bestfit = 1.32;
  // const double ue4_data_bestfit = 0.116;
  // const double um4_data_bestfit = 0.135;

  GridDefinition griddef;
  
  //griddef.define_maya_grid();
  // int mi_stride   = 20; //  5 pts
  // int uei_stride  = 10; //  5 pts
  // int umui_stride = 10; //  5 pts
  
  griddef.define_proper_bounds_course();
  int mi_stride   = 5; //  5 pts
  int uei_stride  = 2; //  5 pts
  int umui_stride = 2; //  5 pts
  
  //const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
  //const int dm2_grdpts(1), ue4_grdpts(25), umu4_grdpts(25);  // for best-fit profiles: dm2 fixed
  //const int dm2_grdpts(25), ue4_grdpts(1), umu4_grdpts(25);  // for best-fit profiles: Ue4 fixed
  //const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(1);  // for best-fit profiles: Um4 fixed
  const int dm2_grdpts( griddef.dm2_grdpts ), ue4_grdpts( griddef.ue4_grdpts ), umu4_grdpts( griddef.umu4_grdpts );
  
  // const double dm2_lowbound(0.01), dm2_hibound(100);
  // //const double dm2_lowbound(dm2_data_bestfit), dm2_hibound(dm2_data_bestfit); // best-fit value
  // //const double ue4_lowbound(0.01), ue4_hibound(0.5);
  // const double ue4_lowbound(ue4_data_bestfit), ue4_hibound(ue4_data_bestfit);   // best-fit value
  // const double umu4_lowbound(0.01), umu4_hibound(0.5);
  //const double umu4_lowbound(um4_data_bestfit), umu4_hibound(um4_data_bestfit);
  const double dm2_lowbound(  griddef.dm2_lowbound ),  dm2_hibound(  griddef.dm2_hibound );
  const double ue4_lowbound(  griddef.ue4_lowbound ),  ue4_hibound(  griddef.ue4_hibound );
  const double umu4_lowbound( griddef.umu4_lowbound ), umu4_hibound( griddef.umu4_hibound );
  
  // tag for some save files
  std::string tag = "DL_full";

  // --------------------------------------------------

  // For profiling
  // we need to align dm2 grid with data-bestfit. we do this by:
  // (1) find the closet grid edge
  // (2) calculate the log(dm^2) shift to align that closest grid to the data bestfit value
  int idx_nearest_dm2 = -1;
  float offset_logdm = 0;
  // if ( dm2_grdpts>1 ) {
  //   float min_offset = 1e9;
  //   float bestfit_logdm = TMath::Log10(sqrt(dm2_data_bestfit));
  //   for (int i=0; i<=dm2_grdpts; i++) {
  //     float grid_logdm = float(i)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound));
  //     float offset = bestfit_logdm-grid_logdm;
  //     if ( fabs(offset)<fabs(min_offset) ) {
  // 	min_offset = offset;
  // 	idx_nearest_dm2 = i;
  //     }
  //   }
  //   offset_logdm = min_offset;
  // }
  std::cout << "nearest dm2 gridpoint: " << idx_nearest_dm2 << std::endl;
  std::cout << "offset logdm: " << offset_logdm << std::endl;  
  
  // --------------------------------------------------
  
  // initialize the minimizer
  SBNllminimizer minimizer( xml );
  std::cout<<"Initialized minimizer"<<std::endl;
  minimizer.use_polar_coords = false; 

  // now get the parameters we are testing
  // specific_entry=0;
  int mi_base   = griddef.binindex_v[specific_entry][0];
  int uei_base  = griddef.binindex_v[specific_entry][1];
  int umui_base = griddef.binindex_v[specific_entry][2];

  int mnu_base = sqrt( griddef.oscpar_v[specific_entry][0] );
  int ue_base  = griddef.oscpar_v[specific_entry][1];
  int um_base  = griddef.oscpar_v[specific_entry][2];
  
  // open output files
  std::ofstream chifile;
  std::string textid = ZeroPadNumber(specific_entry, 5);
  chifile.open("chis_seek_"+textid+".txt", std::ios_base::app);
  
  std::ofstream pseudofile;
  pseudofile.open( "pseudoexp_"+textid+".txt", std::ios_base::app );

  std::ofstream fitresult_file;
  fitresult_file.open( "fitresult_"+textid+".txt", std::ios_base::app );

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

  // following is obsolete now that we rely on Maya's grid
  // // there were 400 mass spectra pregenerated in the original version
  // //  - this makes sure the grid points line up
  // int mi_base_new = mi_base*(400/dm2_grdpts);
  // //mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2
  // float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
  // float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
  // //float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
  // float mnu_base = float(mi_base)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound));
  // if ( dm2_grdpts>1 ) {
  //   mnu_base += offset_logdm;
  // }
  // mnu_base = pow( 10., mnu_base );
  
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
  TrueChi.InitRandomNumberSeeds( (double)specific_entry );
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

  float tot_pseudo = 0.0;
  float tot_gridscan = 0.0;
  float tot_simplex = 0.0;
  float tot_migrad = 0.0;

  // Across several fake experiments for this grid point:
  for(int expi = 0; expi < nFakeExp; expi++){

    // ===================================================
    // START OF PSEUDOEXPERIMENT FIT ROUTINE
    // ===================================================

    auto start_pseudoexp = high_resolution_clock::now();
    
    std::cout<<"EXPERIMENT NUMBER: "<<expi<<std::endl;
    std::vector<float> fakeData = TrueChi.GeneratePseudoExperiment();

    fitresult_file << "----------------------------" << std::endl;    

    // save pseudo-experiment draw to text file
    pseudofile << "[" << expi << "] ";
    fitresult_file << "[" << expi << "] ";
    for ( auto& bincount : fakeData ) {
      pseudofile  << bincount << " ";
      fitresult_file << bincount << " ";
    }
    pseudofile << std::endl;
    fitresult_file << std::endl;

    // option to hardcode fakedata spectrum for testing
    if(false){
      fakeData=SetFakeData(fakeData);
    }

    // get the -2llh at the throw point (pt)
    TMatrixD cov_pt = GetTotalCov(fakeData, oscSpec, *oscFracSys_collapsed);
    float ptfit     = GetLLHFromVector(fakeData, oscSpec, cov_pt, false);
    fitresult_file << "ptfit: " << ptfit << std::endl;

    // ===================================================
    // the fit procedures
    // ===================================================
    
    // option to run traditional grid search
    float chi_min_grid =1000000;

    // now do a grid search to draw confidence levels
    // we want to know the best locations
    class grid_result_t {
    public:
      grid_result_t()
	: llh(0),
	  mi(0),
	  uei(0),
	  umui(0)
      {};
      grid_result_t( float llh_, int mi_, int uei_, int umui_ )
	: llh(llh_),
	  mi(mi_),
	  uei(uei_),
	  umui(umui_)
      {};
      virtual ~grid_result_t() {};
      
      bool operator<( const grid_result_t& rhs ) const
      {
	if ( llh < rhs.llh )
	  return true;
	return false;
      };
      
      float llh;
      int mi;
      int uei;
      int umui;
      
    };

    auto start_gridscan = high_resolution_clock::now();    
    
    std::vector< grid_result_t > result_list_v; /// container to put in results
    int num_sorts = 0;
    
    if(true){

      // strided grid search
      // result is a 5x5x5=125 pt search on the grid
      
      int num_results = (griddef.dm2_grdpts/mi_stride)*(griddef.ue4_grdpts/uei_stride)*(griddef.umu4_grdpts/umui_stride);
      result_list_v.reserve( num_results );

      // float e_app_in = 4*pow(ue_base,2)*pow(um_base,2);
      // float e_dis_in = 4*pow(ue_base,2)*(1-pow(ue_base,2));
      // float m_dis_in = 4*pow(um_base,2)*(1-pow(um_base,2));

      float uei_min, umui_min, dmi_min;

      // loop over all the grid points
      for(int mi_in = 0; mi_in <griddef.dm2_grdpts; mi_in += mi_stride ){
        for(int uei_in = 0; uei_in < griddef.ue4_grdpts; uei_in += uei_stride ){
          for(int umui_in = 0; umui_in < griddef.umu4_grdpts; umui_in += umui_stride ){

	    int kindex = griddef.binindex_to_kindex( mi_in, uei_in, umui_in );
	    
            // int mi_in_new = mi_in*(400/numgridpts);
            // float ue_val = pow(10.,(uei_in/float(numgridpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
            // float um_val = pow(10.,(umui_in/float(numgridpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
            // //float mnu_val = pow(10.,((mi_in+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
	    // float mnu_val = float(mi_in)/float(numgridpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound));
	    // mnu_val += offset_logdm;
	    // mnu_val = pow(10.,mnu_val);

	    float ue_val = griddef.Ue4_bin_v[ uei_in ];
	    float um_val = griddef.Umu4_bin_v[ umui_in ];
	    float mnu_val = sqrt( griddef.dm2_bin_v[ mi_in ] );

            // get the oscillated spectra and covariance matrix
            SBNspec inSpec =  minimizer.getOscSpectra( mnu_val, ue_val, um_val );
            // inSpec.PrintCollapsedVector();
            inSpec.RemoveMCError();
            SBNchi innerChi(inSpec, *covFracSys);
            TMatrixD * inFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
            int c = innerChi.FillCollapsedFractionalMatrix(inFracSys_collapsed);
            TMatrixD cov_grid = GetTotalCov(fakeData, inSpec, *inFracSys_collapsed);
            float chi_test = GetLLHFromVector(fakeData, inSpec, cov_grid, false);
            delete inFracSys_collapsed;

	    // save the result
	    result_list_v.push_back( grid_result_t( chi_test, mi_in, uei_in, umui_in ) );
	    
            if (chi_test<chi_min_grid){
              chi_min_grid = chi_test;
	      dmi_min = mi_in;
	      uei_min = uei_in;
	      umui_min = umui_in;
            }
          }
        }
      }

      std::sort( result_list_v.begin(), result_list_v.end() );
      
      // print the top results
      std::cout << "GRID MINIMUM" << std::endl;
      std::cout<< "chi2=" << chi_min_grid
	       << " binidex=( " << dmi_min << ", " << uei_min << ", " << umui_min << ") "
	       << " par=(" << griddef.dm2_bin_v[dmi_min] << ", " << griddef.Ue4_bin_v[uei_min] << ", " << griddef.Umu4_bin_v[umui_min] << ")"
	       << std::endl;

      std::cout << "GRID BEST POINTS" << std::endl;
      int iresult =0 ;
      for ( auto& result : result_list_v ) {
	std::cout << "[" << iresult << "] chi2=" << result.llh
		  << " grid_index=(" << result.mi << ", " << result.uei << ", " << result.umui << ")"
		  << std::endl;
	iresult++;
      }
      
    }// end of if to run traditional grid search

    auto stop_gridscan = high_resolution_clock::now();
    auto duration_gridscan = duration_cast<milliseconds>(stop_gridscan - start_gridscan);
    
    // use top grid points to seed simplex minimizer
    minimizer.algoName = "Simplex";
    
    class simplex_result_t {
      
    public:
      simplex_result_t()
	: llh(0),
	  mi(0),
	  uei(0),
	  umui(0)      
      {};
      simplex_result_t( float llh_, float mi_, float uei_, float umui_ )
	: llh(llh_),
	  mi(mi_),
	  uei(uei_),
	  umui(umui_)
      {};
      virtual ~simplex_result_t() {};
      
      bool operator<( const simplex_result_t& rhs ) const
      {
	if ( llh < rhs.llh )
	  return true;
	return false;
      };
      
      float llh;
      float mi;
      float uei;
      float umui;
      
    };
    
    std::vector< simplex_result_t > simplex_v;
    std::vector<double> beststart(4,1.0e3);
    double min_minimizer = 1e6;
    
    auto start_simplex = high_resolution_clock::now();
    
    if(DO_SIMPLEX_SEARCH){
      
      std::vector<double> temp;
      
      std::cout<<"starting data fit loop"<<std::endl;
      
      for (int iiresult=0; iiresult<NTOP_GRID_TO_FIT; iiresult++) {
	if ( iiresult>=result_list_v.size() )
	  break;

	auto const& result = result_list_v.at(iiresult);
	
	// int startval = 5*x;
	// int numgridpts=25;
	// float ue_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
	// float um_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
	// float mnu_val = pow(10.,((startval/float(numgridpts))+.5)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound)) );
	int uei_base = result.uei;
	int umui_base = result.umui;
	
	float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
	float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
	float mnu_base = sqrt( griddef.dm2_bin_v[ result.mi ] );

	fitresult_file << result.llh << " " // chi2
		       << result.mi  << " " << result.uei << " " << result.umui << " " // bin indices
		       << griddef.dm2_bin_v[ result.mi ] << " " << griddef.Ue4_bin_v[result.uei] << " " << griddef.Umu4_bin_v[result.umui]; // osc par values
	
	temp = minimizer.doFit( fakeData, mnu_base, ue_base , um_base );

	fitresult_file << " "
		       << temp[0] << " " // post minimizer chi2
		       << temp[1] << " " << temp[2] << " " << temp[3] // fitted osc par values
		       << std::endl;
	
	if (temp[0] < min_minimizer){
	  min_minimizer=temp[0];
	  beststart=temp;
	}
	
	simplex_v.push_back( simplex_result_t( temp[0], sqrt(temp[1]), temp[2], temp[3] ) );
	
      } //end of loop over start points
      
      // std::cout << "PAST MINIMUM"  << std::endl;
      // temp = minimizer.doFit( data_v, sqrt(3.70), 0.6, 0.0 );
      // if ( temp[0] < min_minimizer ) {
      //   min_minimizer = temp[0];
      //   beststart = temp;
      // }

      std::sort( simplex_v.begin(), simplex_v.end() );
      
      std::cout << "SIMPLEX RESULTS" << std::endl;
      int isimplex = 0;
      for (auto& simplex : simplex_v ) {
	std::cout << "[" << isimplex << "] chi2=" << simplex.llh
		  << " mu=" << simplex.mi << " Ue4=" << simplex.uei << " Umu4=" << simplex.umui
		  << std::endl;
      }
      
    }//end of it SIMPLEX SEARCH

    auto stop_simplex = high_resolution_clock::now();
    auto duration_simplex = duration_cast<milliseconds>(stop_simplex - start_simplex);  

    // now run the minimizer
    auto start_migrad = high_resolution_clock::now();  
    minimizer.algoName = "MIGRAD";
    
    std::vector<double> migrad_result
      = minimizer.doFit( fakeData, simplex_v.front().mi, simplex_v.front().uei, simplex_v.front().umui );

    std::cout << "FINAL MIGRAD RESULT" << std::endl;
    std::cout << "[0] chi2=" << migrad_result[0]
	      << " mu=" << migrad_result[1] << " Ue4=" << migrad_result[2] << " Umu4=" << migrad_result[3]
	      << std::endl;
    fitresult_file << migrad_result[0] << " "
		   << migrad_result[1] << " "
		   << migrad_result[2] << " "
		   << migrad_result[3]
		   << std::endl;
    
    auto stop_migrad = high_resolution_clock::now();
    auto duration_migrad = duration_cast<milliseconds>(stop_migrad - start_migrad);  
    
    // SAVE BEST FIT
    if ( simplex_v.front().llh < migrad_result[0] ) {
      // save simplex output instead
      beststart[0] = simplex_v.front().llh;
      beststart[1] = simplex_v.front().mi*simplex_v.front().mi;
      beststart[2] = simplex_v.front().uei;
      beststart[3] = simplex_v.front().umui;
    }
    else {
      beststart = migrad_result;
    }

    auto stop_pseudoexp = high_resolution_clock::now();
    auto duration_pseudoexp = duration_cast<milliseconds>(stop_pseudoexp - start_pseudoexp);

    fitresult_file << "time: "
		   << duration_pseudoexp.count() << " "
		   << duration_gridscan.count() << " "
		   << duration_simplex.count() << " "            
		   << duration_migrad.count()
		   << std::endl;

    tot_pseudo += duration_pseudoexp.count();
    tot_gridscan += duration_gridscan.count();
    tot_simplex += duration_simplex.count();
    tot_migrad += duration_migrad.count();
    
    // // trying a variety of starts
    // double min_minimizer =100000000;
    // std::vector<float> beststart;
    // if(true){
    //   std::cout<<"starting loop"<<std::endl;
    //   double temp;
    //   // std::vector<float> m_v={0.01,2,8};
    //   // std::vector<float> u_v={0.01,0.08,0.4};
    //   for(int x=0;x<5;x++){
    //         int startval = 5*x;
    //         int numgridpts=25;
    //         float ue_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
    //         float um_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
    //         float mnu_val = pow(10.,((startval*(500/float(numgridpts))+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));

    //         temp = minimizer.doFit( fakeData, mnu_val,ue_val,um_val)[0];
    //         if (temp < min_minimizer){
    //           min_minimizer=temp;
    //           // beststart={m_v[x],u_v[y],u_v[z]};
    //       //   }
    //       // }
    //     }
    //   } //end of loop over start points
    // }    
    chifile << expi << " " << ptfit << " " << beststart[0] << " " << beststart[1] << " " << beststart[2] << " " << beststart[3] << std::endl;
    
  } //end of fake experiment loop

  fitresult_file << "---------------------------------" << std::endl;
  fitresult_file << "total times: "
		 << tot_pseudo << " "
		 << tot_gridscan << " "
		 << tot_simplex << " "            
		 << tot_migrad
		 << std::endl;
  
  
  std::cout << "DONE!!" << std::endl;
  fitresult_file.close();
  
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
