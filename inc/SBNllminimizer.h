#ifndef _SBN_LL_MINIMIZER_H__
#define _SBN_LL_MINIMIZER_H__

#include "SBNgenerate.h"
#include "SBNspec.h"
#include "prob.h"

namespace sbn {

  /**
   * @brief Osc likelihood minimizer
   *
   * Uses:
   *  - a cache-enabled SBNgenerate instance to produce oscillated spectrum
   *  - Katie's negative log-likelihood calculation
   *  - ROOT minimizer framework: https://root.cern.ch/root/html534/tutorials/fit/NumericalMinimization.C.html
   */
  class SBNllminimizer {

  public:

    SBNllminimizer( std::string );
    virtual ~SBNllminimizer();

    NeutrinoModel _osc_model; // stores the oscillation parameters
    SBNgenerate _gen;     // makes predicted spectra
    SBNchi* _chi;         // makes covariance matrix, calculates chi-2, generates pseudo-experiments
    int nBins_e; // want to not hard code this info
    int nBins_mu; // want to not hard code this info
    int nBins; // want to not hard code this info
    bool use_polar_coords; // fit in polar form for U: r=|Ue4|^2+|Um4|^2 and phi=atan2(|U_e4|,|Um4|)
    std::string algoName; // options: simplex, migrad, seek

    static double negative_likelihood_ratio( const double* par );
    std::vector<double> doFit( std::vector<float>& obs_bins, float dm_start,float ue_start, float um_start );

    TMatrixD* covFracSys;


    static TMatrixD GetTotalCov(const std::vector<float>& obsSpec,
				const SBNspec& expSpec,
				const TMatrixD& Mfracsys);
    static float GetLLHFromVector(const std::vector<float>& obsSpec,
				  const SBNspec& expSpec,
				  const TMatrixD& Msys,
				  bool prints);
    static void UtoPolar( const float Ue4, const float Um4, float& r, float& phi );
    static void PolarToU( const float r, const float phi, float& Ue4, float& Um4 );

    std::vector<float> _observed_bins;
    void setObservedBinValues( std::vector<float>& obs ) { _observed_bins = obs; };

    SBNspec getOscSpectra( float dm41, float ue4, float um4 );

  protected:

    int _niters; // steps we've taken
    static SBNllminimizer* _active_copy;

  };

  typedef  double (SBNllminimizer::*SBNllminizerFCN)( const double* par );
}

#endif
