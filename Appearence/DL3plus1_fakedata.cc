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

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

// define helper functions
void printbinedges();
SBNspec GetOscillatedSpectra(SBNspec cvSpec, SBNspec massSpec,
			float e_app, float e_dis, float m_dis);
TMatrixD GetTotalCov(const std::vector<float>& obsSpec,const SBNspec& expSpec,const TMatrixD& Mfracsys);
float GetLLHFromVector(const std::vector<float>& obsSpec, const SBNspec& expSpec,const TMatrixD& Msys,bool prints);
void testfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
std::string ZeroPadNumber(int num, int digits);

// define some global variables
std::string xml = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/xml/TotalThreePlusOne_full.xml";
bool gen = false;
bool printbins = true;
int mass_start = -1;
std::string tag = "DL_full";
// set these parameters at the very start
const double dm2_lowbound(0.01), dm2_hibound(100);
const double ue4_lowbound(0.01), ue4_hibound(0.5);
const double umu4_lowbound(0.01), umu4_hibound(0.5);
// to genergate dm2_grdpts = 100
const int dm2_grdpts(25), ue4_grdpts(25), umu4_grdpts(25);
const int nBins_e(12),nBins_mu(19);
const int nBins = nBins_e+nBins_mu;
TMatrixD cov_osc(nBins,nBins);
TMatrixD cov_cv(nBins,nBins);
double  mnu, ue, umu;
int count;
SBNspec  appSpec, innerSpec, oscSpec;
SBNspec cvSpec("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/"+tag+".SBNspec.root",xml);
std::array<double,nBins>  a_pred;
std::ofstream coordfile;
std::ofstream chifile;
std::ofstream chi_sinee_file;
std::ofstream chi_sinmu_file;
std::ofstream chi_sinemu_file;
std::ofstream specfile;
Float_t z[5],x[5],y[5],errorz[5];
std::vector<std::tuple<SBNspec,double>> a_sinsqSpec;
TMatrixD * covFracSys_collapsed;
std::vector<float> data_v;
TFile * fsys = new TFile("/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot.SBNcovar.root","read");
TMatrixD * covFracSys = (TMatrixD*)fsys->Get("frac_covariance");
bool sin2profile = false;
bool sin2slice = false;

int main(int argc, char* argv[]){
	// get the input integer: require 1
	int fakeset = atoi(argv[1]);
	int sintype = atoi(argv[2]);
	std::string textid = ZeroPadNumber(fakeset, 1);

	// open output text files
	coordfile.open("bins_fakedata.txt", std::ios_base::app);
	chifile.open("chis_fakedata"+textid+".txt", std::ios_base::app);
	specfile.open("spec_fakedata2.txt", std::ios_base::app);
	// Print binedges for easier plotting
	if(printbins) printbinedges();

	// Load up the necesary bits and store them in vectors on the stack **
	cvSpec.Scale("fullosc",0.0);
	cvSpec.CollapseVector();
	cvSpec.RemoveMCError();
	specfile<<"cvspec:"<<std::endl;
	for(int i=0;i<nBins;i++){
		specfile<<cvSpec.collapsed_vector[i]<<" ";
	}
	specfile<<std::endl;
	// cvSpec.PrintFullVector();
	// make a tuple of the spectra and mass term
	for(int mi = 0; mi < 400; mi++){
		mnu = pow(10.,((mi+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
		std::stringstream stream;
		stream << std::fixed << std::setprecision(4) << 2*log10(mnu);
		std::string infile = "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/MassSpectraBigBins/"+tag+"_SINSQ_dm_"+stream.str()+".SBNspec.root";
		auto inspec = SBNspec(infile,xml);
		// inspec.Scale("ext",0.0);	// since we're subtracting this spectrum, we want to make sure we're not subtracting the background.
		inspec.CollapseVector();
		inspec.RemoveMCError();
		std::tuple<SBNspec,double> singletup (inspec,mnu);
		a_sinsqSpec.push_back(singletup);
	}

  // manually setting the data vector
  // set1
  if(fakeset==1) data_v={20,12,6,6,4,5,0,1,0,1,2,1,18,128,201,279,236,312,328,314,229,194,177,155,144,148,121,89,63,37,48};
  // set2
  if(fakeset==2) data_v={4,5,3,2,5,5,1,1,1,3,7,4,23,185,292,334,353,400,462,435,372,305,279,228,190,190,133,119,75,60,54};

  // now do a grid search to draw confidence levels
	int idx=0;
  float ue_min, um_min, m_min;
  float grid_min = 99999;
	for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
		for(int uei_base = 0; uei_base < ue4_grdpts; uei_base++){
			for(int umui_base = 0; umui_base < umu4_grdpts; umui_base++){
				std::cout<<idx<<std::endl;
				idx+=1;
					//mnu =m41, ue4 = sqrt(.5sin^2(2theta14)), um4 = sqrt(.5sin^2(2theta24)), sin2 term = 4(ue4**2)(um4**2)
				int mi_base_new = mi_base*(400/dm2_grdpts);
				float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
				float um_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
				float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				// calculate scaling factors
				float e_app = 4*pow(ue_base,2)*pow(um_base,2);
				float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
				float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));
				// current test model
				// std::cout << "NU Base Model: m41^2:" << mnu_base << " ue:" << ue_base << " um:" << um_base << std::endl;
				// std::cout <<"scale factors: "<<e_app<<" "<<e_dis<<" "<<m_dis<<std::endl;

				// get the oscillated spectra and cov matrix
				oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
				SBNchi TrueChi(cvSpec, *covFracSys);
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

				chifile<<llh_osc<<std::endl;
			}// end of loop over base umu
		}//end of loop over base ue
	}//end of loop over base mass

  std::cout<< grid_min <<" "<<ue_min <<" "<<um_min<<" "<<m_min ;


  // first get the best fit using a minimizer
  // will be the first line of the chis file
  Int_t ierflg = 0;
  Double_t dm2_step = .01;
  TString dm2_name = "m41";
  Double_t ue4_step = 0.001;
  TString ue4_name = "ue4";
  Double_t umu4_step = 0.001;
  TString umu4_name = "umu4";
  Double_t arglist[2];
  Double_t dm2_val,dm2_err,dm_low,dm_high;
  Int_t dm2_num;
  Double_t ue4_val,ue4_err,ue4_low,ue4_high;
  Int_t ue4_num;
  Double_t umu4_val,umu4_err,umu4_low,umu4_high;
  Int_t umu4_num;
  Int_t npari,nparx,istat;
  // Now ready for minimization step 0 = max steps, 1=tolerance
  // give a large maximum number of steps (tends to take ~125 calls)
  arglist[0] = 500;
  arglist[1] = 3;

  // next initialize the minimizer
  TMinuit * gMinuit_simple= new TMinuit(3);
  gMinuit_simple->SetFCN(testfcn);
  Double_t dm2_start = m_min;
  Double_t ue4_start =ue_min;
  Double_t umu4_start = um_min;
  // start values near lower bound, but some wiggle room
  // initialize the parameters
  gMinuit_simple->mnparm(0,dm2_name,dm2_start,dm2_step,sqrt(dm2_lowbound),sqrt(dm2_hibound),ierflg);
  gMinuit_simple->mnparm(1,ue4_name,ue4_start,ue4_step,ue4_lowbound,ue4_hibound,ierflg);
  gMinuit_simple->mnparm(2,umu4_name,umu4_start,umu4_step,umu4_lowbound,umu4_hibound,ierflg);
  // run the minimizer
  gMinuit_simple->mnexcm("SEEk", arglist ,2,ierflg);
  // get the parameter outputs
  Double_t fmin,fedm,errdef;
  gMinuit_simple->mnpout(0,dm2_name,dm2_val,dm2_err,dm_low,dm_high,dm2_num);
  gMinuit_simple->mnpout(1,ue4_name,ue4_val,ue4_err,ue4_low,ue4_high,ue4_num);
  gMinuit_simple->mnpout(2,umu4_name,umu4_val,umu4_err,umu4_low,umu4_high,umu4_num);
  gMinuit_simple->mnstat(fmin,fedm,errdef,npari,nparx,istat);
  // std::cout<<"pt, fmin, coarse grid, fine grid"<<std::endl;
  std::cout<<fmin<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<std::endl;
  chifile<<fmin<<" "<<dm2_val<<" "<<ue4_val<<" "<<umu4_val<<std::endl;

	// float ue4_val = 0;
	// float umu4_val = .21;
	// float dm2_val = sqrt(3.37);
	// finally get some stuff to print...
	float e_app_best = 4*pow(ue4_val,2)*pow(umu4_val,2);
	float e_dis_best = 4*pow(ue4_val,2)*(1-pow(ue4_val,2));
	float m_dis_best = 4*pow(umu4_val,2)*(1-pow(umu4_val,2));
	std::cout<<e_app_best<<" "<<e_dis_best<<" "<<m_dis_best<<std::endl;
	if(sin2profile){
		chi_sinee_file.open("chis_sinee_fakedata"+textid+".txt", std::ios_base::app);
		chi_sinmu_file.open("chis_sinmu_fakedata"+textid+".txt", std::ios_base::app);
		chi_sinemu_file.open("chis_sinemu_fakedata"+textid+".txt", std::ios_base::app);
		// perform fit over sin2 terms insteds
		if(sintype==0){
			for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
				// start with sin2ee
				for(int sinei_base = 0; sinei_base < 25; sinei_base++){
					float llh_best=9999999;
					std::cout<<mi_base<<" "<<sinei_base<<std::endl;
					for(int sinmui_base = 0; sinmui_base < 25; sinmui_base++){
						for(int sinemui_base = 0; sinemui_base < 25; sinemui_base++){
							int mi_base_new = mi_base*(400/dm2_grdpts);
							float sine_base = pow(10.,(sinei_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float sinmu_base = pow(10.,(sinmui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float sinemu_base = pow(10.,(sinemui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
							// calculate scaling factors
							float e_app = sinemu_base;
							float e_dis = sine_base;
							float m_dis = sinmu_base;

							// get the oscillated spectra and cov matrix
							oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
							SBNchi TrueChi(cvSpec, *covFracSys);
							TrueChi.ReloadCoreSpectrum(&oscSpec);
							TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
							int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
							cov_osc= GetTotalCov(data_v, oscSpec, *oscFracSys_collapsed);
							// arguments: obs,exp,cov
							float llh_osc = GetLLHFromVector(data_v, oscSpec, cov_osc,false);
							if(llh_osc<llh_best){
								llh_best=llh_osc;
							}
						}// end of loop over sin_emu
					}//end of loop over sin_mu
					chi_sinee_file<<llh_best-0<<std::endl;
				}//end of loop over sin_e
			}// end of loop over mass
		}
		// perform fit over sin2 terms insteds
		if(sintype==1){
			for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
				// next with sin2mu
				for(int sinmui_base = 0; sinmui_base < 25; sinmui_base++){
					float llh_best=9999999;
					std::cout<<mi_base<<" "<<sinmui_base<<std::endl;
					for(int sinei_base = 0; sinei_base < 25; sinei_base++){
						for(int sinemui_base = 0; sinemui_base <25; sinemui_base++){
							int mi_base_new = mi_base*(400/dm2_grdpts);
							float sine_base = pow(10.,(sinei_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float sinmu_base = pow(10.,(sinmui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float sinemu_base = pow(10.,(sinemui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
							// calculate scaling factors
							float e_app = sinemu_base;
							float e_dis = sine_base;
							float m_dis = sinmu_base;

							// get the oscillated spectra and cov matrix
							oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
							SBNchi TrueChi(cvSpec, *covFracSys);
							TrueChi.ReloadCoreSpectrum(&oscSpec);
							TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
							int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
							cov_osc= GetTotalCov(data_v, oscSpec, *oscFracSys_collapsed);
							// arguments: obs,exp,cov
							float llh_osc = GetLLHFromVector(data_v, oscSpec, cov_osc,false);
							if(llh_osc<llh_best){
								llh_best=llh_osc;
							}
						}// end of loop over sin_emu
					}//end of loop over sin_mu
					chi_sinmu_file<<llh_best-0<<std::endl;
				}//end of loop over sin_e
			}// end of loop over mass
		}


		// perform fit over sin2 terms insteds
		if(sintype==2){
			for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
				// finally sin2emu
				for(int sinemui_base = 0; sinemui_base < 25; sinemui_base++){
					float llh_best=9999999;
					std::cout<<mi_base<<" "<<sinemui_base<<std::endl;
					for(int sinmui_base = 0; sinmui_base < 25; sinmui_base++){
						for(int sinei_base = 0; sinei_base <25; sinei_base++){
							int mi_base_new = mi_base*(400/dm2_grdpts);
							float sine_base = pow(10.,(sinei_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float sinmu_base = pow(10.,(sinmui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float sinemu_base = pow(10.,(sinemui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
							float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
							// calculate scaling factors
							float e_app = sinemu_base;
							float e_dis = sine_base;
							float m_dis = sinmu_base;

							// get the oscillated spectra and cov matrix
							oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
							SBNchi TrueChi(cvSpec, *covFracSys);
							TrueChi.ReloadCoreSpectrum(&oscSpec);
							TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
							int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
							cov_osc= GetTotalCov(data_v, oscSpec, *oscFracSys_collapsed);
							// arguments: obs,exp,cov
							float llh_osc = GetLLHFromVector(data_v, oscSpec, cov_osc,false);
							if(llh_osc<llh_best){
								llh_best=llh_osc;
							}
						}// end of loop over sin_emu
					}//end of loop over sin_mu
					chi_sinemu_file<<llh_best-0<<std::endl;
				}//end of loop over sin_e
			}// end of loop over mass
		}

	}

	if(sin2slice){
		chi_sinee_file.open("chis_sinee_slice_fakedata"+textid+".txt", std::ios_base::app);
		chi_sinmu_file.open("chis_sinmu_slice_fakedata"+textid+".txt", std::ios_base::app);
		chi_sinemu_file.open("chis_sinemu_slice_fakedata"+textid+".txt", std::ios_base::app);
		// perform fit over sin2 terms insteds
		for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
			// start with sin2ee
			for(int sinei_base = 0; sinei_base < 25; sinei_base++){
				float llh_best=9999999;
				std::cout<<mi_base<<" "<<sinei_base<<std::endl;

				int mi_base_new = mi_base*(400/dm2_grdpts);
				float sine_base = pow(10.,(sinei_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
				float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				// calculate scaling factors
				float e_app=0;
				float e_dis = sine_base;
				float m_dis =0;
				if(fakeset==1){
					e_app = 0.217;
					m_dis = 0.746;
				}
				if(fakeset==2){
					e_app = 0.00003;
					m_dis = 0.237;
				}
				// get the oscillated spectra and cov matrix
				oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
				SBNchi TrueChi(cvSpec, *covFracSys);
				TrueChi.ReloadCoreSpectrum(&oscSpec);
				TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
				int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
				cov_osc= GetTotalCov(data_v, oscSpec, *oscFracSys_collapsed);
				// arguments: obs,exp,cov
				float llh_osc = GetLLHFromVector(data_v, oscSpec, cov_osc,false);

				chi_sinee_file<<llh_osc-0<<std::endl;
			}//end of loop over sin_e
		}// end of loop over mass

		// perform fit over sin2 terms insteds

		for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
			// next with sin2mu
			for(int sinmui_base = 0; sinmui_base < 25; sinmui_base++){
				float llh_best=9999999;
				std::cout<<mi_base<<" "<<sinmui_base<<std::endl;

				int mi_base_new = mi_base*(400/dm2_grdpts);
				float sinmu_base = pow(10.,(sinmui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
				float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				// calculate scaling factors
				float e_app = 0;
				float e_dis = 0;
				float m_dis = sinmu_base;
				if(fakeset==1){
					e_app = 0.217;
					e_dis = 0.684;
				}
				if(fakeset==2){
					e_app = 0.00003;
					e_dis = 0.00048;
				}

				// get the oscillated spectra and cov matrix
				oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
				SBNchi TrueChi(cvSpec, *covFracSys);
				TrueChi.ReloadCoreSpectrum(&oscSpec);
				TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
				int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
				cov_osc= GetTotalCov(data_v, oscSpec, *oscFracSys_collapsed);
				// arguments: obs,exp,cov
				float llh_osc = GetLLHFromVector(data_v, oscSpec, cov_osc,false);

				chi_sinmu_file<<llh_osc<<std::endl;
			}//end of loop over sin_e
		}// end of loop over mass

		// perform fit over sin2 terms insteds
		for(int mi_base = 0; mi_base <dm2_grdpts; mi_base++){
			// finally sin2emu
			for(int sinemui_base = 0; sinemui_base < 25; sinemui_base++){
				float llh_best=9999999;
				std::cout<<mi_base<<" "<<sinemui_base<<std::endl;
				int mi_base_new = mi_base*(400/dm2_grdpts);
				float sinemu_base = pow(10.,(sinemui_base/float(25)*TMath::Log10(1/.0001) + TMath::Log10(0.0001)));
				float mnu_base = pow(10.,((mi_base_new+.5)/float(400)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
				// calculate scaling factors
				float e_app = sinemu_base;
				float e_dis = 0;
				float m_dis = 0;
				if(fakeset==1){
					m_dis = 0.746;
					e_dis = 0.684;
				}
				if(fakeset==2){
					m_dis = 0.237;
					e_dis = 0.00048;
				}

				// get the oscillated spectra and cov matrix
				oscSpec =  GetOscillatedSpectra(cvSpec, std::get<0>(a_sinsqSpec.at(mi_base_new)), e_app, e_dis, m_dis);
				SBNchi TrueChi(cvSpec, *covFracSys);
				TrueChi.ReloadCoreSpectrum(&oscSpec);
				TMatrixD * oscFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
				int run = TrueChi.FillCollapsedFractionalMatrix(oscFracSys_collapsed);
				cov_osc= GetTotalCov(data_v, oscSpec, *oscFracSys_collapsed);
				// arguments: obs,exp,cov
				float llh_osc = GetLLHFromVector(data_v, oscSpec, cov_osc,false);

				chi_sinemu_file<<llh_osc<<std::endl;
			}//end of loop over sin_e
		}// end of loop over mass


	}


	NeutrinoModel testModel(dm2_val, sqrt(.5), sqrt(.5));
	SBNgenerate * gen = new SBNgenerate(xml,testModel);
	gen->WritePrecomputedOscSpecs("fake1");
	std::stringstream stream;
	stream << std::fixed << std::setprecision(4) << 2*log10(dm2_val);
	std::string infile = "fake1_SINSQ_dm_"+stream.str()+".SBNspec.root";
	auto inspec = SBNspec(infile,xml);
	SBNspec bestSpec =  GetOscillatedSpectra(cvSpec, inspec, e_app_best, e_dis_best, m_dis_best);
	SBNchi bestChi(bestSpec, *covFracSys);
	bestChi.ReloadCoreSpectrum(&bestSpec);
	TMatrixD * bestFracSys_collapsed =(TMatrixD*)fsys->Get("frac_covariance");
	int run_best = bestChi.FillCollapsedFractionalMatrix(bestFracSys_collapsed);

	TMatrixD bestFracSys_tot= GetTotalCov(data_v, bestSpec, *bestFracSys_collapsed);
	TMatrixD invcov = bestFracSys_tot;
	invcov.Invert();
	float chisqTest_e = 0;
	for(int i = 0; i < nBins_e; i++){
		for(int j = 0; j < nBins_e; j++){
			chisqTest_e += (data_v[i] - bestSpec.collapsed_vector[i])*invcov[i][j]*(data_v[j] -bestSpec.collapsed_vector[j]);
		}
	}
	std::cout<<"chi2test_e: "<<chisqTest_e<<std::endl;

	float chisqTest_m = 0;
	for(int i = nBins_e; i < nBins; i++){
		for(int j = nBins_e; j < nBins; j++){
			chisqTest_m += (data_v[i] - bestSpec.collapsed_vector[i])*invcov[i][j]*(data_v[j] -bestSpec.collapsed_vector[j]);
		}
	}
	std::cout<<"chi2test_m: "<<chisqTest_m<<std::endl;


	std::cout<<"best spec:"<<std::endl;
	bestSpec.PrintFullVector();
	std::cout<<std::endl;
	bestSpec.PrintCollapsedVector();

	specfile<<"best cov:"<<std::endl;
	for(int i=0;i<nBins;i++){
		specfile<<(*bestFracSys_collapsed)(i,i)<<" ";
	}
	specfile<<std::endl;

	return 0;
} // end of main function


//-----------------------------------------------------------------------------
//----------------------HELPER FUNCTIONS---------------------------------------
//----------------------------------------------------------------------------

void printbinedges(){
	// funtion that prints the bins  to the output textfile
	for(int mi = 0; mi < dm2_grdpts; mi++){
		mnu = pow(10.,((mi+.5)/float(dm2_grdpts)*TMath::Log10(sqrt(dm2_hibound)/sqrt(dm2_lowbound)) + TMath::Log10(sqrt(dm2_lowbound))));
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
	// std::cout<<"chi2: "<<chisqTest<<" det: "<<log(tempcov.Determinant())<<std::endl;
	chisqTest += log(Msys.Determinant());

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
