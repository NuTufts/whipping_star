#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
// #include "TMinuitMinimizer.h"
#include "TMinuit.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "SBNcovariance.h"
#include "prob.h"
#include "assert.h"
#include "SBNllminimizer.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

// define helper functions
void printbinedges();
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
std::string ZeroPadNumber(int num, int digits);

// define some global variables
std::string xml = "/cluster/tufts/wongjiradlabnu/kmason03/mergewfermi/whipping_star/xml/TotalThreePlusOne_full.xml";
bool gen = false;
bool printbins = true;
int mass_start = -1;
std::string tag = "DL_full";
// set these parameters at the very start

// uboone best fit
// const double dm2_data_bestfit = 3.69939;
// const double ue4_data_bestfit = 0.5;
// const double um4_data_bestfit = 0.01;

// global best fit
const double dm2_data_bestfit = 1.32;
const double ue4_data_bestfit = 0.116;
const double um4_data_bestfit = 0.135;

//const double dm2_lowbound_fixed(0.01), dm2_hibound_fixed(100);

// // config for compatibility with Maya's grid
//const double dm2_lowbound(0.0101158), dm2_hibound(112.202);
const double dm2_lowbound(0.01), dm2_hibound(100.0);
const double ue4_lowbound(0.01), ue4_hibound(0.5);
const double umu4_lowbound(0.01), umu4_hibound(0.5);
const int dm2_grdpts(100), ue4_grdpts(50), umu4_grdpts(50); // profile with dm2 fixed

// // config for profile with dm2 fixed
// const double dm2_lowbound(dm2_data_bestfit), dm2_hibound(dm2_data_bestfit); // profile with dm2 fixed
// const double ue4_lowbound(0.01), ue4_hibound(0.5);
// const double umu4_lowbound(0.01), umu4_hibound(0.5);
// const int dm2_grdpts(1), ue4_grdpts(25), umu4_grdpts(25); // profile with dm2 fixed

// // config for profile with Ue4 fixed
// const double dm2_lowbound(0.01), dm2_hibound(100.0);
// const double ue4_lowbound(0.5), ue4_hibound(0.5);
// const double umu4_lowbound(0.01), umu4_hibound(0.5);
// const int dm2_grdpts(25), ue4_grdpts(1), umu4_grdpts(25); // profile with Ue4 fixed

// config for profile with Um4 fixed
// const double dm2_lowbound(0.01), dm2_hibound(100.0);
// const double ue4_lowbound(0.01), ue4_hibound(0.5);
// const double umu4_lowbound(um4_data_bestfit), umu4_hibound(um4_data_bestfit);
// const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(1); // profile with Ue4 fixed

// const double ue4_lowbound(0.01), ue4_hibound(0.5);
// const double umu4_lowbound(0.01), umu4_hibound(0.5);
// // to genergate dm2_grdpts = 100
// const int dm2_grdpts(1), ue4_grdpts(25), umu4_grdpts(25); // profile with dm2 fixed
// //const int dm2_grdpts(25), ue4_grdpts(1), umu4_grdpts(25); // profile with Ue4 fixed

const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;
TMatrixD cov_osc(nBins,nBins);
TMatrixD cov_cv(nBins,nBins);
double  mnu, ue, umu;
int count;
SBNspec  appSpec, innerSpec, oscSpec;
SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/DL_full_Bkg.SBNspec.root",xml);
std::array<double,nBins>  a_pred;
std::ofstream coordfile;
std::ofstream chifile;
std::ofstream chi_sinee_file;
std::ofstream chi_sinmu_file;
std::ofstream chi_sinemu_file;
std::ofstream specfile;
std::ofstream gridptfile;
Float_t z[5],x[5],y[5],errorz[5];
std::vector<std::tuple<SBNspec,double>> a_sinsqSpec;
TMatrixD * covFracSys_collapsed;
std::vector<float> data_v;
TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
bool sin2slice = false;
std::string pseudo_experiment_file = "/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/data/pseudo_experiments.csv";
std::vector<float> get_pseudo_experiment_data( int idnum, std::string inputfile_path,std::vector<float>& oscpars ); // forward declaration

int main(int argc, char* argv[]){

  int pseudo_exp_idnum = -1;

  if ( argc==2 ) {
    pseudo_exp_idnum = std::atoi( argv[1] );
    std::cout << "//////////////////////////////////////////" << std::endl;
    std::cout << "PSEUDO-EXPERIMENT ID SPECIFIED: " << pseudo_exp_idnum << std::endl;
    std::cout << "//////////////////////////////////////////" << std::endl;    
  }

  // open output text files
  coordfile.open("bins_data.txt", std::ios_base::app);
  chifile.open("chis_data.txt", std::ios_base::app);
  specfile.open("spec_data.txt", std::ios_base::app);
  gridptfile.open("gridpts_data.txt", std::ios_base::app);
  if(printbins) printbinedges();

  // stop after printing bin edges
  if (true) {
    coordfile.close();    
    return 0;
  }
  
  // Load up the necesary bits and store them in vectors on the stack **
  cvSpec.Scale("fullosc",0.0);
  cvSpec.CollapseVector();
  cvSpec.RemoveMCError();
  specfile<<"cvspec:"<<std::endl;
  for(int i=0;i<nBins;i++){
    specfile<<cvSpec.collapsed_vector[i]<<" ";
  }
  specfile<<std::endl;

  // manually setting the data vector
  std::vector<float> file_oscpar_v;  
  if ( pseudo_exp_idnum==-1 ) {
    // default measured data
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "SET DATA TO KATIE'S DLLEE SELECTION=" << pseudo_exp_idnum << " ]" << std::endl;    
    data_v={4., 1., 1., 2., 5., 3., 8., 0., 1., 0., 4., 2., 26., 192., 276. ,401., 389., 463., 439., 482., 395., 353., 303., 233., 240., 177., 118., 109., 109., 85., 58.};
    for ( auto const& bindata : data_v )
      std::cout << bindata << " ";
    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;    
  }
  else {
    data_v = get_pseudo_experiment_data( pseudo_exp_idnum, pseudo_experiment_file, file_oscpar_v );
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "PSEUDO-EXPERIMENT DATA[ ID=" << pseudo_exp_idnum << " ]" << std::endl;
    for ( auto const& bindata : data_v )
      std::cout << bindata << " ";
    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
  }

  // set data to NULL expectation, for sensitivity curves
  //data_v={1.77973,3.299,2.8411,2.84015,3.06677,2.90407,2.99971,2.52731,2.45532,2.19775,
  //	  6.17792,2.0456,30.4971,193.128,279.326,376.987,339.661,413.51,439.355,406.634,
  //	  367.536,327.103,307.313,262.759,217.072,188.681,148.804,114.063,88.3917,71.2635,64.0484}; // 31 bins


  // load in the full systematic matrix - created on fermilab grid
  TFile * fsys = new TFile( "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
  TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");

  // create sbnchi - used to get cov matrix at given osc + throw universes
  SBNchi cvChi(cvSpec, *covFracSys);
  // save cv covar for plotting purposes
  TMatrixD * cvFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
  int a = cvChi.FillCollapsedFractionalMatrix(cvFracSys_collapsed);

  // initialize the minimizer
  // we'll use it to find the best fit and also to generate spectra at different points on the grid
  SBNllminimizer minimizer( xml );
  std::cout<<"Initialized minimizer"<<std::endl;

  // we need to align dm2 grid with data-bestfit. we do this by:
  // (1) find the closet grid edge
  // (2) calculate the log(dm^2) shift to align that closest grid to the data bestfit value
  int idx_nearest_dm2 = -1;
  float offset_logdm = 0;
  if ( dm2_grdpts>1 && false) {
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
  
  // now do a grid search to draw confidence levels
  int idx=0;
  float ue_min, um_min, m_min;
  float grid_min = 99999;
  for(int mi_base = 0; mi_base < dm2_grdpts; mi_base++){
    for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
      for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){
	std::cout<<idx<<std::endl;

	float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
	float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
	//float mnu_base = ((mi_base)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound)));
	float mnu_base = ((mi_base+0.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound)));
	if ( dm2_grdpts>1 ) {
	  mnu_base += offset_logdm;
	}
	mnu_base = pow(10.,mnu_base);
	
	// calculate scaling factors
	float e_app = 4*pow(ue_base,2)*pow(um_base,2);
	float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
	float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));

	// get the oscillated spectra as the expectation for this grid point
	SBNspec oscSpec =  minimizer.getOscSpectra( mnu_base, ue_base, um_base );
	oscSpec.RemoveMCError();
	
	// calculate chi2
	SBNchi TrueChi(oscSpec, *covFracSys);
	TrueChi.ReloadCoreSpectrum(&oscSpec);
	TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
	int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);

	// arguments: obs,exp,cov
	cov_osc= GetTotalCov(data_v, oscSpec, *oscFracSys_collapsed);
	// arguments: obs,exp,cov
	float llh_osc = GetLLHFromVector(data_v, oscSpec, cov_osc,false);
	if(llh_osc<grid_min){
	  grid_min = llh_osc;
	  ue_min = ue_base;
	  um_min = um_base;
	  m_min = mnu_base;
	}

	gridptfile << idx << " " << ue_base << " " << um_base << " " << mnu_base << " " << llh_osc << std::endl;
	
	chifile<<llh_osc<<std::endl;
	idx += 1;
      }// end of loop over base umu
    }//end of loop over base ue
  }//end of loop over base mass

  std::cout<< grid_min <<" "<<ue_min <<" "<<um_min<<" "<<m_min ;


  // trying a variety of starts
  double min_minimizer =100000000;
  std::vector<double> beststart;
  if(true){
    std::cout<<"starting loop"<<std::endl;
    std::vector<double> temp;
    for(int x=0;x<5;x++){
      int startval = 5*x;
      int numgridpts=25;
      float ue_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
      float um_val = pow(10.,(startval/float(numgridpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
      float mnu_val = pow(10.,((startval*(500/float(numgridpts))+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
      temp = minimizer.doFit( data_v, mnu_val,ue_val,um_val);
      if (temp[0] < min_minimizer){
	min_minimizer=temp[0];
	beststart=temp;
	//   }
	// }
      }
    } //end of loop over start points
  }
  chifile<<beststart[0]<<" "<<beststart[1]<<" "<<beststart[2]<<" "<<beststart[3]<<std::endl;

  // float ue4_val = 0;
  // float umu4_val = .21;
  // float dm2_val = sqrt(3.37);
  // finally get some stuff to print...
  float e_app_best = 4*pow(beststart[2],2)*pow(beststart[3],2);
  float e_dis_best = 4*pow(beststart[2],2)*(1-pow(beststart[2],2));
  float m_dis_best = 4*pow(beststart[3],2)*(1-pow(beststart[3],2));
  std::cout<<e_app_best<<" "<<e_dis_best<<" "<<m_dis_best<<std::endl;
  return 0;
} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

void printbinedges(){
  // funtion that prints the bins  to the output textfile
  for(int mi = 0; mi < dm2_grdpts; mi++){
    //mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
    mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) ));
    //mnu = pow(10.,((mi)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
    coordfile << pow(mnu,2) << " ";
  }
  coordfile << std::endl;
  
  for(int uei = 0; uei < ue4_grdpts; uei++){
    ue = pow(10.,(uei/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
    coordfile << ue << " ";
  }
  coordfile << std::endl;
  
  for(int umui = 0; umui < umu4_grdpts; umui++){
    umu = pow(10.,(umui/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
    coordfile << umu << " ";
  }
  coordfile << std::endl;
  return;
}//end of print bins function


SBNspec GetOscillatedSpectra(SBNspec cvSpec, SBNspec massSpec,
			     float e_app, float e_dis, float m_dis){
  // function to take the cv spec and return the oscillated spec
  // this is only set up to do all osc at once.
  // inputs:
  // cvspec: cv spectum
  // massspec: oscillated spec for the given mass, maximum oscillations
  // e_app,e_dis,m_dis oscillation scaling parameters
  // output oscspec: oscillated spectrum
  SBNspec testSpec = cvSpec;
  testSpec.Scale("fullosc",0.0);
  massSpec.Scale("fullosc",e_app);
  massSpec.Scale("bnb",-1*m_dis);
  massSpec.Scale("nue",-1*e_dis);
  massSpec.Scale("ext",0.0);
  // massSpec.PrintCollapsedVector();
  testSpec.Add(&massSpec);
  // std::cout<<"oscillated spec"<<std::endl;
  // testSpec.PrintCollapsedVector();
  return testSpec;
}//end of GetOscillatedSpectra

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
      // add in stat errors start with CNP for "data" errors
      if(i==j){
	if (expSpec.collapsed_vector[i] >0 ){
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
	float chisqTest;

	// inv cov for chi2calc
	TMatrixD invcov = Msys;
	invcov.Invert();

	// add the chi2-like part
	chisqTest = 0;
	for(int i = 0; i < nBins; i++){
		for(int j = 0; j < nBins; j++){
			// (obsi-predi)*(invcov)*(obsj-predj)
			if(i==j && prints) std::cout<<i<<" "<<obsSpec[i]<<" "<<expSpec.collapsed_vector[i]<<" "<<Msys[i][j]<<" "<<((obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]))<<std::endl;
			chisqTest += (obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]);
		}
	}
	// now need ln(det(2Pi*M))
	// TMatrixD tempcov = 2*3.14159265358979323846*Msys;
	// std::cout<<"chi2: "<<chisqTest<<" det: "<<log(Msys.Determinant())<<std::endl;
	chisqTest += log(Msys.Determinant());

  // if (chisqTest<207){
  //   std::cout<<"chi2: "<<chisqTest-log(Msys.Determinant())<<" det: "<<log(Msys.Determinant())<<std::endl;
  //   std::cout<<"best cov:"<<std::endl;
  //   for(int i=0;i<nBins;i++){
  //     std::cout<<(Msys[i][i])<<" ";
  //   }
  //   std::cout<<std::endl;
  //   assert (1==2);
  // }

	return chisqTest;
}//end of GetLLHFromSpectra

void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

	// function to calculate the chi2 between fake universe and the mc events with osc weights
	float e_app = 4*pow(par[1],2)*pow(par[2],2);
	float e_dis = 4*pow(par[1],2)*(1-pow(par[1],2));
	float m_dis = 4*pow(par[2],2)*(1-pow(par[2],2));

	// find the closest mass spectra
	int lowidx, highidx;
	double prevval = 0;
	for(int i = 1;i<a_sinsqSpec.size();i++){
		// std::cout<<std::get<1>(a_sinsqSpec.at(i))<<" "<<par[0]<<" "<<prevval<<std::endl;
		if (par[0]==std::get<1>(a_sinsqSpec.at(i))){
			lowidx = i;
			highidx =i;
			break;
		}
		else if (par[0]<std::get<1>(a_sinsqSpec.at(i)) && par[0] >prevval ){
			lowidx = i-1;
			highidx =i;
			break;
		}
		else if( i == a_sinsqSpec.size()-1){
			lowidx = i;
			highidx =i;
		}
		else prevval = std::get<1>(a_sinsqSpec.at(i));
	}

	int closeidx = lowidx;
	if (lowidx <0) lowidx =0;
	if (lowidx<a_sinsqSpec.size()-2){
		double diff1 = par[0] - std::get<1>(a_sinsqSpec.at(lowidx));
		double diff2 = std::get<1>(a_sinsqSpec.at(highidx)) - par[0];
		if (diff2 <diff1) closeidx =highidx;
	}
	// get spectra and covar
	SBNspec inspec = std::get<0>(a_sinsqSpec.at(closeidx));
	SBNspec newSpec = GetOscillatedSpectra(cvSpec, inspec,e_app,e_dis, m_dis);
	SBNchi tmpChi(newSpec, *covFracSys);
	TMatrixD * tmpFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
	int b = tmpChi.FillCollapsedFractionalMatrix(tmpFracSys_collapsed);
	TMatrixD tmpcov = GetTotalCov(data_v, newSpec, *tmpFracSys_collapsed);
	// calculate -2LLH
	f =  GetLLHFromVector(data_v, newSpec, tmpcov, false);
	delete tmpFracSys_collapsed;
} // end of fcn function

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

std::vector<float> get_pseudo_experiment_data( int idnum, std::string inputfile_path,
					       std::vector<float>& oscpars )
{
  std::cout << "//// get_pseudo_experiment_data idnum=" << idnum << " ///////" << std::endl;
  
  std::ifstream infile( inputfile_path.c_str() );
  std::string s;
  while ( !infile.eof() ) {
    s = "";
    std::getline(infile,s);
    if ( s=="" )
      break;
    
    // get first number which is pseudo-experiment index
    int lineid = -1;
    std::string::size_type pos1 = 0;
    std::string::size_type pos2 = 0;    
    pos2 = s.find(',',pos1);
    if ( pos2!=std::string::npos ) {
      std::string sub = s.substr(pos1,pos2-pos1);
      try {
	lineid = std::atoi( sub.c_str() );
      }
      catch (std::exception& e) {
	continue;	
      }
    }

    if ( lineid!=idnum )
      continue;

    std::cout << "found idnum match: " << lineid << std::endl;

    // if it is right idnum, continue to parse line
    int tokensread = 0;
    int nbins = 0;
    int ibin = 0;
    pos1 = pos2+1;
    pos2 = s.find(",",pos1);
    while ( pos2!=std::string::npos ) {
      std::string sub = s.substr(pos1,pos2-pos1);
      if ( tokensread==0 ) {
	// first token (after id num) is the number of bins
	nbins = std::atoi( sub.c_str() );
	data_v.resize(nbins,0);
      }
      else if ( tokensread>0 && ibin<nbins ) {
	float binval = std::atof( sub.c_str() );
	data_v[ibin] = binval;
	ibin++;
      }
      else if ( ibin==nbins ) {
	float parval = std::atof( sub.c_str() );
	oscpars.push_back(parval);
      }
      tokensread++;
      pos1 = pos2+1;
      pos2 = s.find(",",pos1);
      std::cout << "read token  pos1=" << pos1 << " pos2=" << pos2 << std::endl;
    }
    std::cout << "parsed line" << std::endl;
    break;
  }// while loop

  return data_v;
}
