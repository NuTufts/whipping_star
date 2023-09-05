#include "SBNllminimizer.h"

#include <string>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

namespace sbn {

  SBNllminimizer* SBNllminimizer::_active_copy = NULL;

  std::vector<double> SBNllminimizer::llh_components_v = { 0, 0, 0};

  SBNllminimizer::SBNllminimizer( std::string xml_config )
    : _osc_model(0.01,0.707,0.707), // these values are arbitrary
      _gen( xml_config, _osc_model, true ),
      _chi( NULL ),
      nBins_e(12), // this needs to be not hard-coded
      nBins_mu(19),
      nBins(31),
      use_polar_coords(false),
      fit_only_one_par(false), // if true, we fix the parameters for all but one parameter
      unfixed_par_index(-1),
      fix_one_par(false),
      fixed_par_index(-1),
      algoName("Simplex"), //Migrad, Seek      
      covFracSys(NULL)
  {
    // set as active copy if the first instance made
    if ( !_active_copy )
      _active_copy = this;
    
    // get the config obj through SBNgen
    TiXmlDocument doc(_gen.xmlname.c_str());
    bool loadOkay = doc.LoadFile();
    std::string otag = "SBNllminimizer::SBNllminimizer\t||\t";
    if(!loadOkay){
      std::cout<<otag<<"ERROR: Failed to load XML configuration file: "<<_gen.xmlname<<std::endl;
      std::cout<<otag<<"ERROR: This generally means broken .xml brackets or attribute syntax."<<std::endl;
      exit(EXIT_FAILURE);
    }
    TiXmlHandle hDoc(&doc);
    // look at preloaded list
    TiXmlElement *ppreload = doc.FirstChildElement("Precalced");
    if ( !ppreload ) {
      std::cout << otag << "Filed to load required precalced covariance matrix." << std::endl;
      std::cout << otag << "Need it listed in a <Precalced> block in the xml file" << std::endl;
      std::cout << otag << "looked in: " << _gen.xmlname << std::endl;
      exit(EXIT_FAILURE);
    }
    else {
      TiXmlElement *pquantity = ppreload->FirstChildElement("quantity");
      while(pquantity){
	std::string quantity_type = pquantity->Attribute("type");
	if ( quantity_type=="covar" ) {
	  std::string cov_fpath = pquantity->Attribute("filename");
	  std::string cov_name  = pquantity->Attribute("name");
	  TFile _fsys(cov_fpath.c_str(),"read");
	  covFracSys = (TMatrixD*)_fsys.Get(cov_name.c_str()); //"frac_covariance"
	  _fsys.Close();
	  // check the matrix, it should be NxN where N=_gen.num_bins_total
	  if ( covFracSys->GetNcols()!=_gen.num_bins_total || covFracSys->GetNrows()!=_gen.num_bins_total ) {
	    std::cout << otag << "loaded total covar matrix does not have the expected shape" << std::endl;
	    std::cout << otag << "  cov.shape=(" << covFracSys->GetNcols() << "," << covFracSys->GetNrows() << ")" << std::endl;
	    std::cout << otag << "  num_bins_total=" << _gen.num_bins_total << " (all subchannel bins, before collapsed)" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  else {
	    std::cout << otag << "Total (all subchannel) covar loaded. "
		      << "cov.shape=(" << covFracSys->GetNcols() << "," << covFracSys->GetNrows() << ")"
		      << std::endl;
	  }
	}
	pquantity = pquantity->NextSiblingElement("quantity");
      }
    }//end of preloaded xml element

    _chi = new SBNchi( _gen.spec_central_value, covFracSys );
    _chi->is_verbose = false;

  }

  SBNllminimizer::~SBNllminimizer()
  {
    if ( covFracSys )
      delete covFracSys;
    if ( _chi )
      delete _chi;
  }

  void SBNllminimizer::UtoPolar( const float Ue4, const float Um4, float& r, float& phi )
  {
    r = TMath::Sqrt(Ue4*Ue4 + Um4*Um4);
    phi = TMath::ATan2( Um4, Ue4 );
  }

  void SBNllminimizer::PolarToU( const float r, const float phi, float& Ue4, float& Um4 )
  {
    Ue4 = r*TMath::Cos(phi);
    Um4 = r*TMath::Sin(phi);
  }
  
  // static function
  double SBNllminimizer::negative_likelihood_ratio( const double* par )
  {
    //std::cout << _active_copy << std::endl;

    // function to calculate the chi2 between fake universe and the mc events with osc weights
    float logdm2 = par[0];
    float dm     = sqrt( exp(logdm2) ); // gross

    float Ue4    = par[1];
    float Um4    = par[2];
    if ( _active_copy->use_polar_coords ) {
      // convert from polar coordinates
      SBNllminimizer::PolarToU( par[1], par[2], Ue4, Um4 );
    }

    float e_app  = 4*pow(Ue4,2)*pow(Um4,2);     // sin^2(2theta_mue)
    float e_dis  = 4*pow(Ue4,2)*(1-pow(Ue4,2)); // sin^2(2theta_ee)
    float m_dis  = 4*pow(Um4,2)*(1-pow(Um4,2)); // sin^2(2theta_mumu)

    //std::cout << "dm=" << dm << " logdm^2=" << logdm2 << " Ue4=" << Ue4 << " Um4=" << Um4 << std::endl;

    // update the osc model
    _active_copy->_osc_model = NeutrinoModel( dm, Ue4, Um4 );
    _active_copy->_gen.regenerate_osc( _active_copy->_osc_model ); ///< regenerates spectrum with new dm4x

    // now we scale the different components
    _active_copy->_gen.spec_osc_sinsq.Scale("fullosc",e_app);
    _active_copy->_gen.spec_osc_sinsq.Scale("bnb",-1*m_dis);
    _active_copy->_gen.spec_osc_sinsq.Scale("nue",-1*e_dis);
    _active_copy->_gen.spec_osc_sinsq.Scale("ext",0.0);
    _active_copy->_gen.spec_central_value.Scale("fullosc",0.0);

    // std::cout << "pre-osc full central value -------- " << std::endl;
    // _active_copy->_gen.spec_central_value.PrintFullVector(true);
    // std::cout << "pre-scale oscillated value -------- " << std::endl;
    // _active_copy->_gen.spec_osc_sinsq.PrintFullVector(true);
    // std::cout << "----------------------------------- " << std::endl;

    _active_copy->_gen.spec_central_value.Add(&(_active_copy->_gen.spec_osc_sinsq));
    //_active_copy->_gen.spec_central_value.PrintFullVector(true);

    // pass oscillated spectrum prediction to SBNchi instance
    // this will store bin values and also build covariance matrix for this prediction
    _active_copy->_gen.spec_central_value.RemoveMCError();
    _active_copy->_chi->ReloadCoreSpectrum( &(_active_copy->_gen.spec_central_value) );

    // get fractional inverse cov matrix
    //TFile * fsys = new TFile( "/cluster/tufts/wongjiradlabnu/kmason03/whipping_star/data/systematics/katieversion_bigbins_tot_copy.SBNcovar.root","read");
    //TMatrixD * inv_frac_cov =(TMatrixD*)fsys->Get("frac_covariance");
    //SBNchi oscChi( _active_copy->_gen.spec_central_value, *inv_frac_cov);
    SBNchi oscChi( _active_copy->_gen.spec_central_value, *(_active_copy->covFracSys) ); // pass-by-copy AND copied in constructor
    TMatrixD inv_frac_cov;
    int b = oscChi.FillCollapsedFractionalMatrix( &inv_frac_cov ); // modifies and fills the matrix with collapsed frac sys

    // this is varying from call to call as the relative contributions are modified by oscillations
    // changing the flux more or less
    // std::cout<<"temp cov: "<<logdm2*logdm2 << " " << Ue4 << " " << Um4 << std::endl;
    // for(short i = 0; i < 31; i++){
    //   for(short j = 0; j < 31; j++){
    //     std::cout<< inv_frac_cov(i,j)<<" ";
    //   }
    //   std::cout<<std::endl;
    // }

    // scalce frac inverse cov matrix to expectation and include CNP stat error to the diag
    TMatrixD tmpcov = SBNllminimizer::GetTotalCov(_active_copy->_observed_bins,
						  _active_copy->_gen.spec_central_value,
						  inv_frac_cov );
    

    // calculate -2LLH
    double f = SBNllminimizer::GetLLHFromVector(_active_copy->_observed_bins,
						_active_copy->_gen.spec_central_value,
						tmpcov,
						false);
    
    if ( std::isnan(f) || std::isinf(f) ) {
      std::cout << "SBNllminimizer::negative_likelihood_ratio NLLR is bad = " << f << std::endl;
      throw std::runtime_error("bad likelihood calculated");
    }
    if ( _active_copy->_niters%100==0 )
      std::cout << "SBNllminimizer::negative_likelihood_ratio iter=" << _active_copy->_niters << " NLLR=" << f << std::endl;
    _active_copy->_niters++;

    //delete fsys;
    //delete inv_frac_cov;
    return f;
  }

  std::vector<double> SBNllminimizer::doFit( std::vector<float>& obs_bins, float dm_start,float ue_start, float um_start )
  {
    std::string minName =  "Minuit";

    _active_copy = this;
    _active_copy->setObservedBinValues( obs_bins ); // we need a number of bins check

    ROOT::Math::Minimizer* min =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);


    ROOT::Math::Functor f(&SBNllminimizer::negative_likelihood_ratio,3);
    min->SetFunction(f);

    const float dm2_lowbound(1.0e-3), dm2_hibound(1e3);
    // const float ue4_lowbound(0.0),  ue4_hibound(0.70710678118);
    // const float umu4_lowbound(0.0), umu4_hibound(0.70710678118);
    const float ue4_lowbound(0.0),  ue4_hibound(0.5);
    const float umu4_lowbound(0.0), umu4_hibound(0.5);
    float logdm2_low  = log(dm2_lowbound);
    float logdm2_high = log(dm2_hibound);

    // min->SetVariable( 0, "log(dm^2)", (dm_start), 1 );
    std::cout << "Initial log_dm2 step: " << fabs(logdm2_high-logdm2_low)/3.0 << std::endl;
    min->SetVariable( 0, "log(dm^2)", log(dm_start*dm_start), fabs(logdm2_high-logdm2_low)/3.0 );
    min->SetVariableLimits( 0, logdm2_low, logdm2_high );
    
    if ( !use_polar_coords ) {
      // standard Ue4 and Um4 space
      min->SetVariable( 1, "Ue4", ue_start, 0.1 );
      min->SetVariable( 2, "Um4", um_start, 0.1 );
      min->SetVariableLimits( 1, ue4_lowbound, ue4_hibound );
      min->SetVariableLimits( 2, umu4_lowbound, umu4_hibound );

      if ( fit_only_one_par ){
	for (int p=0; p<3; p++) {
	  if ( p!=unfixed_par_index ) {
	    min->SetVariableValue(p,fixed_par_values[p]);
	    min->FixVariable(p);
	  }
	  else {
	    min->ReleaseVariable(p);
	  }
	}
      }
      else if ( fix_one_par && fixed_par_index>=0 & fixed_par_index<3) {
	min->SetVariableValue( fixed_par_index, fixed_par_values[fixed_par_index] );
	min->FixVariable( fixed_par_values[fixed_par_index] );
      }
    }
    else {
      // polar (Ue4,Um4) space to help enforce bounds: |Ue4|^2+|Um4|^2<=1/2
      float rstart = TMath::Sqrt(ue_start*ue_start+um_start*um_start);
      float phistart = atan2(um_start,ue_start);
      min->SetVariable( 1, "r", rstart, 0.1 );
      min->SetVariable( 2, "phi", phistart, 0.1 );
      min->SetVariableLimits( 1, 0, TMath::Sqrt(0.5) );
      min->SetVariableLimits( 2, 0, 1.57079632679 ); // 0 to pi/2
    }
    
    min->SetMaxFunctionCalls(1000); // for Minuit/Minuit2
    min->SetMaxIterations(1000);  // for GSL
    min->SetTolerance(0.01);
    min->SetPrintLevel(1);

    _active_copy->_niters = 0;
    min->Minimize();
    const double* results = min->X();
    std::vector<double> bestfit(4,0.0);
    bestfit[0] = min->MinValue(); // chi2
    bestfit[1] = exp(results[0]); // dm^2
    if ( !use_polar_coords ) {
      bestfit[2] = results[1];
      bestfit[3] = results[2];
    }
    else {
      float ue4min = 0.;
      float um4min = 0.;
      SBNllminimizer::PolarToU( results[1], results[2], ue4min, um4min );
      bestfit[2] = ue4min;
      bestfit[3] = um4min;
    }

    std::cout << "MINIMIZER RESULTS:" << std::endl;
    std::cout << "  Chi^2 = " << bestfit[0] << std::endl;
    std::cout << "  dm^2 = " << bestfit[1] << " eV^2" << std::endl;
    std::cout << "  dm   = " << sqrt(bestfit[1]) << " eV" << std::endl;
    std::cout << "  Ue4 = "  << bestfit[2] << std::endl;
    std::cout << "  Um4 = "  << bestfit[3] << std::endl;
    
    // destroy minimizer
    delete min;
    
    /// return the results
    return bestfit;
  }

  TMatrixD SBNllminimizer::GetTotalCov(const std::vector<float>& obsSpec, const SBNspec& expSpec, const TMatrixD& Mfracsys) {
    // function to take the fractional Msys and return totl Msys+Mstat
    // inputs:
    // obsSpec: "data" spectra
    // predSpec: "MC" spectra
    // Mfracsys: fractional (flux+xsec+detvar) covariance matrix
    int nBins = _active_copy->nBins;
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


  float SBNllminimizer::GetLLHFromVector(const std::vector<float>& obsSpec,
					 const SBNspec& expSpec,
					 const TMatrixD& Msys,
					 bool prints)
  {
    // // function to calculate a chi2 (shape + rate)
    // // inputs:
    // // obsSpec: "data" vector
    // // predSpec: "MC" spectra
    // // Mfracsys: total (flux+xsec+detvar) covariance matrix
    double chisqTest = 0;
    int nBins = _active_copy->nBins;
    // expSpec->RemoveMCError();

    // inv cov for chi2calc
    TMatrixD invcov = Msys;
    invcov.Invert();

    llh_components_v.resize(3,0.0);

    // add the chi2-like part
    chisqTest = 0;
    for(int i = 0; i < nBins; i++){
      for(int j = 0; j < nBins; j++){
	// (obsi-predi)*(invcov)*(obsj-predj)
    	if(i==j && false) std::cout<<i<<" "
				   <<obsSpec[i]<<" "
				   <<expSpec.collapsed_vector[i]<<" "
				   <<Msys[i][j]<<" "
				   <<((obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]))
				   <<std::endl;
    	chisqTest += (obsSpec[i] - expSpec.collapsed_vector[i])*invcov[i][j]*(obsSpec[j] - expSpec.collapsed_vector[j]);
      }
    }
    llh_components_v[1] = chisqTest;
    // now need ln(det(2Pi*M))
    // TMatrixD tempcov = 2*3.14159265358979323846*Msys;
    // std::cout<<"chi2: "<<chisqTest<<" det: "<<log(tempcov.Determinant())<<std::endl;
    double MDet = log(Msys.Determinant());
    llh_components_v[2] = MDet;
    chisqTest += MDet;

    llh_components_v[0] = chisqTest;

    return chisqTest;
  }//end of GetLLHFromSpectra

  SBNspec SBNllminimizer::getOscSpectra( float dm4x, float Ue4, float Um4 )
  {
    _osc_model = NeutrinoModel( dm4x, Ue4, Um4 );
    _gen.regenerate_osc( _osc_model ); ///< regenerates spectrum with new dm4x

    float e_app = 4*pow(Ue4,2)*pow(Um4,2);     // sin^2(2theta_mue)
    float e_dis = 4*pow(Ue4,2)*(1-pow(Ue4,2)); // sin^2(2theta_ee)
    float m_dis = 4*pow(Um4,2)*(1-pow(Um4,2)); // sin^2(2theta_mumu)

    // now we scale the different components
    _gen.spec_osc_sinsq.Scale("fullosc",e_app);
    _gen.spec_osc_sinsq.Scale("bnb",-1*m_dis);
    _gen.spec_osc_sinsq.Scale("nue",-1*e_dis);
    _gen.spec_osc_sinsq.Scale("ext",0.0);
    _gen.spec_central_value.Scale("fullosc",0.0);
    _gen.spec_central_value.Add(&(_gen.spec_osc_sinsq));

    return _gen.spec_central_value;
  }

}
