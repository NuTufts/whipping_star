#ifndef __KATIE_OPT_PROCEDURE_H__
#define __KATIE_OPT_PROCEDURE_H__
#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>
#include <string>

#include "TRandom3.h"
#include "TMatrixD.h"

#include "SBNllminimizer.h"

#include "GridDefinition.h"

namespace oscfit {

  // has to be header only
  class KatieOptProcedure {
  public:

    KatieOptProcedure()
      : DO_GRID_CHI2_CALC(true),
	DO_SIMPLEX_SEARCH(true),
	DO_DATA_FIT(true),
	SKIP_FULL_RES_GRIDFIT(false),
	_rand(nullptr)
    {
      _rand = new TRandom3(314159);
    };
    virtual ~KatieOptProcedure() {
      delete _rand;
      _rand = nullptr;
    };

    bool DO_GRID_CHI2_CALC;
    bool DO_SIMPLEX_SEARCH;
    bool DO_DATA_FIT;
    bool SKIP_FULL_RES_GRIDFIT;

    // struct for recording grid NLL calculations
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

    // struct holding simplex results
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

    TRandom3* _rand;

    std::vector< std::vector<int> > get_grid_sample_points( int reduce_dim,
							    std::string reduce_method,
							    std::vector<int> fixed_gridpt,
							    GridDefinition& griddef )
    {

      int mi_stride   = 20; //  5 pts
      int uei_stride  = 10; //  5 pts
      int umui_stride = 10; //  5 pts

      std::vector< std::vector<int> > index_list;

      if ( reduce_dim<0 ) {
	// not removing any dims
	for(int mi_base = 0; mi_base < griddef.dm2_grdpts; mi_base += 1 ) {
	  if ( mi_base%mi_stride!=0 )
	    continue;
	  
	  for(int uei_base = 0; uei_base < griddef.ue4_grdpts; uei_base += 1 ) {
	    if ( uei_base%uei_stride!=0 )
	      continue;
	    
	    for(int umui_base = 0; umui_base < griddef.umu4_grdpts; umui_base += 1 ) {
	      if ( umui_base%umui_stride!=0 )
		continue;

	      std::vector<int> index = { mi_base, uei_base, umui_base };
	      index_list.push_back(index);
	    }
	  }
	}
      }
      else {
	// grid scan when we reduce dimensions
	if ( reduce_dim==0 ) {
	  for(int mi_base = 0; mi_base < griddef.dm2_grdpts; mi_base += 1 ) {
	    if ( mi_base%mi_stride!=0 )
	      continue;
	    std::vector<int> index = { mi_base, fixed_gridpt[1], fixed_gridpt[2] };
	    index_list.push_back(index);
	  }
	}
	else if ( reduce_dim==1 ) {
	  for(int uei_base = 0; uei_base < griddef.ue4_grdpts; uei_base += 1 ) {
	    if ( uei_base%uei_stride!=0 )
	      continue;
	    std::vector<int> index = { fixed_gridpt[0], uei_base, fixed_gridpt[2] };
	    index_list.push_back(index);
	  }
	}
	else if ( reduce_dim==2 ) {
	  for(int umui_base = 0; umui_base < griddef.umu4_grdpts; umui_base += 1 ) {
	    if ( umui_base%umui_stride!=0 )
	      continue;	  
	    std::vector<int> index = { fixed_gridpt[0], fixed_gridpt[1], umui_base };
	    index_list.push_back(index);
	  }
	}	    	  
      }
      
      return index_list;
      
    }

    std::vector<double> getReducedLikelihood( std::vector<float>& data_v, std::vector<double>& oscpar,
					      GridDefinition& griddef, sbn::SBNllminimizer& minimizer,
					      int reduce_dim=-1, std::string reduce_method="slice",
					      float slice_par_value=-1,
					      float num_marginal_samples=1000 )
    {
      
      double mnu_base = oscpar[0]; // wants dm
      double ue_base  = oscpar[1]; // wants Ue4
      double um_base  = oscpar[2]; // wants Um4

      minimizer.setObservedBinValues( data_v );      
      
      std::vector<double> NLL_v(3,0);
      
      if ( reduce_method=="profile") {
	// we are removing a par from the likelihood by profiling
	
	// we set the minimizer to only fit one par with simplex. good luck ...
	double par_grid[3] = { log(mnu_base*mnu_base), ue_base, um_base };
	
	minimizer.onlyFitOnePar( reduce_dim, par_grid );
	minimizer.algoName = "Simplex";
	
	// we get the likelihood at this grid point by fitting
	std::vector<double> fit_results = minimizer.doFit( data_v, mnu_base, ue_base , um_base );
	
	// run one more time in case best fit pars not used for last call
	double par_again[3] = { log(fit_results[1]), fit_results[2], fit_results[3] };
	double nll = sbn::SBNllminimizer::negative_likelihood_ratio( par_again );
	NLL_v = sbn::SBNllminimizer::llh_components_v;
      }
      else if ( reduce_method=="marginalize" ) {
	// we throw values of the parameter and average ...
	// minimizer NLL function takes pars such that
	// [0] log(dm2)
	// [1] Ue4
	// [2] Um4
	double par_grid[3] = { log(mnu_base*mnu_base), ue_base, um_base };
	
	for (int v=0; v<3; v++)
	  NLL_v[v] = 0.;
	
	std::cout << "Sampling p(theta) to perform marginalization" << std::endl;

	// sample from prior
	std::vector<double> rand_par(num_marginal_samples,0.0);
	if ( reduce_dim==0 ) {
	  // this should be log-normal? yes if our prior is on scale
	  for (int n=0; n<num_marginal_samples; n++)
	    rand_par[n] = _rand->Uniform(log(griddef.dm2_lowbound),log(griddef.dm2_hibound));
	}
	else if (reduce_dim==1 ) {
	  for (int n=0; n<num_marginal_samples; n++)
	    rand_par[1] = _rand->Uniform( griddef.ue4_lowbound, griddef.ue4_hibound);
	}
	else if ( reduce_dim==2 ) {
	  for (int n=0; n<num_marginal_samples; n++)
	    rand_par[2] = _rand->Uniform( griddef.umu4_lowbound, griddef.umu4_hibound);
	}

	// MC integration by sampling
	for (int n=0; n<num_marginal_samples; n++) {
	  
	  // // sample from prior
	  // if ( reduce_dim==0 ) {
	  //   // this should be log-normal? yes if our prior is on scale
	  //   par_grid[0] = _rand->Uniform(log(griddef.dm2_lowbound),log(griddef.dm2_hibound));
	  // }
	  // else if (reduce_dim==1 ) {
	  //   par_grid[1] = _rand->Uniform( griddef.ue4_lowbound, griddef.ue4_hibound);
	  // }
	  // else if ( reduce_dim==2 ) {
	  //   par_grid[2] = _rand->Uniform( griddef.umu4_lowbound, griddef.umu4_hibound);
	  // }

	  // set osc par
	  par_grid[reduce_dim] = rand_par[n];
	  // calculate NLL
	  float sample_NLL = sbn::SBNllminimizer::negative_likelihood_ratio( par_grid );
	  for (int v=0; v<3; v++)
	    NLL_v[v] += sbn::SBNllminimizer::llh_components_v[v];

	  if ( n%100==0 ) {
	    if ( reduce_dim>0 )
	      std::cout << " [" << n << "] " << par_grid[reduce_dim] << " ~ p(theta) / NLL=" << sample_NLL << std::endl;
	    else
	      std::cout << " [" << n << "] " << exp(par_grid[reduce_dim]) << " ~ p(dm2) / NLL=" << sample_NLL << std::endl;
	  }
	}
	
	// divide by the number of samples to get final NLL
	for (int v=0; v<3; v++)
	  NLL_v[v] /= float(num_marginal_samples);
	std::cout << "Final Marginalized NLL: " << NLL_v[0] << std::endl;
      }
      else if ( reduce_method=="slice" ) {
	// we slice through the grid
      }
    
    
      return NLL_v;
    };


    
    std::vector<double> analyze( std::vector<float>& data_v, std::vector<int>& grid_indices,
				 sbn::SBNllminimizer& minimizer, GridDefinition& griddef,
				 int reduce_dim=-1, std::string reduce_method="slice",
				 float slice_par_value=-1,
				 float num_marginal_samples=1000) {
      /**
	 This produces R for all points of our osc par grid.
	 This requires calculating the NLL at all grid points AND finding the best-fit parameters.

	 if reduce_dim>=0 [must be 0,1,2], we remove that dimension of the par grid defined in griddef.
	 [0]: dm2
	 [1]: Ue4
	 [2]: Um4

	 Methods to reduce the dimension: 
	 'slice': find the nearest grid point nearest to slice_par_value. must be provided.
	 'profile': we minimize the NLL only for the parameter to be profiled
	 'marginalize': we draw N=num_marginal_samples values for the parameter and take the average NLL value
      */

      int mi_stride   = 20; //  5 pts
      int uei_stride  = 10; //  5 pts
      int umui_stride = 10; //  5 pts  
      
      float ue_min = 0;
      float um_min = 0;
      float m_min = 0;
      float grid_min = 99999;

      int NUM_TOP_FOR_SIMPLEX = 5;

      // open log files
      // std::ofstream chifile;
      // chifile.open("chis_data.txt", std::ios_base::app);
      
      // std::ofstream gridptfile;
      // gridptfile.open("gridpts_data.txt", std::ios_base::app);

      // anticipate the number of results on the grid to hold
      int max_results = int(0.01*griddef.dm2_grdpts*griddef.ue4_grdpts*griddef.umu4_grdpts);
      if (max_results<10 )
	max_results = 10;

      minimizer.setObservedBinValues( data_v );

      int idx = 0;
      int num_sorts = 0; // track for perfomance debugging
      std::vector< grid_result_t > result_list_v; /// container to put in results
      result_list_v.reserve( max_results );
      
      auto start_gridscan = std::chrono::high_resolution_clock::now();


      // we sweep part of the grid to find places to seed the fitter
      std::vector< std::vector<int> > gridpt_v = get_grid_sample_points( reduce_dim, reduce_method, grid_indices, griddef );
  
      for ( auto const& indices : gridpt_v ) {

	int mi_base   = indices[0];
	int uei_base  = indices[1];
	int umui_base = indices[2];

	//std::cout<<idx<<std::endl;
	// translate grid point into values needed for optimizer and sbnfit
	float mnu_base = sqrt( griddef.dm2_bin_v[ mi_base ] );	  
	float ue_base  = griddef.Ue4_bin_v[ uei_base ];
	float um_base  = griddef.Umu4_bin_v[ umui_base ];
	
	//std::cout << "grid eval: " << mi_base << " " << uei_base << " " << umui_base << std::endl;
	      
	// calculate the osc scaling factors
	float e_app = 4*pow(ue_base,2)*pow(um_base,2);
	float e_dis = 4*pow(ue_base,2)*(1-pow(ue_base,2));
	float m_dis = 4*pow(um_base,2)*(1-pow(um_base,2));


	// now we need the likelihood
	// we are trying to do two things at once.
	//  (1) we are trying to sample the grid to seed the fitter
	//  (2) we might also be recording the likelihood (reduced or not) on the grid
	
	// calculate the NLL for seeding the fitter
	double par_grid[3] = { log(mnu_base*mnu_base), ue_base, um_base };
	double NLL_seed    = sbn::SBNllminimizer::negative_likelihood_ratio( par_grid );
	std::vector<double> nll_components = sbn::SBNllminimizer::llh_components_v;			      
	
	double NLL = 0;
		
	// we are recoding likelihood values and not just seeding the fitter
	if ( reduce_dim>=0 && !SKIP_FULL_RES_GRIDFIT ) {
	  // record a reduced LL
	  std::vector<double> inputpar = { mnu_base, ue_base, um_base };
	  std::vector<double> reduce_NLL = getReducedLikelihood( data_v, inputpar,
								 griddef, minimizer,
								 reduce_dim, reduce_method,
								 slice_par_value,
								 num_marginal_samples );
	  NLL = reduce_NLL[0];
	  for (int v=0; v<3; )
	    nll_components[v] = reduce_NLL[v];
	}	      
	else {
	  // no dimension reduction or using slice method. just pass the seed NLL
	  NLL = NLL_seed;
	}
	      
	      
	if(NLL<grid_min){
	  grid_min = NLL;
	  ue_min = ue_base;
	  um_min = um_base;
	  m_min = mnu_base;
	}
	
	// log value
	// gridptfile << idx << " "
	// 		 << mi_base << " " << uei_base << " " << umui_base << " "
	// 		 << ue_base << " " << um_base << " " << mnu_base*mnu_base << " "
	// 		 << NLL << " "
	// 		 << nll_components[0] << " " << nll_components[1] << " " << nll_components[2]
	// 		 << std::endl;
	
	
	// on the subgrid, so register grid_result_t into result_list_v
	if ( result_list_v.size() < max_results ) {
	  // result vector not full
	  result_list_v.push_back( grid_result_t( NLL_seed, mi_base, uei_base, umui_base ) );
	  if ( result_list_v.size()==max_results ) {
	    std::sort( result_list_v.begin(), result_list_v.end() );
	    num_sorts++;
	  }
	}
	else {
	  // result vector at capacity
	  if ( result_list_v.back().llh > NLL_seed ) {
	    // add this and sort
	    result_list_v.pop_back(); // remove bad element
	    result_list_v.push_back( grid_result_t( NLL_seed, mi_base, uei_base, umui_base ) );
	    std::sort( result_list_v.begin(), result_list_v.end() );
	    num_sorts++;
	  }
	}

	// log chi2
	//chifile << NLL << " " << nll_components[0] << " " << nll_components[1] << std::endl;
	idx++;
	
      }//end of gridpt list
      
      std::cout << "GRID MINIMUM" << std::endl;
      std::cout<< grid_min <<" "<<ue_min <<" "<<um_min<<" "<<m_min << std::endl;

      if ( num_sorts==0 ) {
	std::sort( result_list_v.begin(), result_list_v.end() );
	num_sorts++;
      }
      
      std::cout << "GRID BEST POINTS" << std::endl;
      int iresult =0 ;
      for ( auto& result : result_list_v ) {
	std::cout << "[" << iresult << "] chi2=" << result.llh
		  << " grid_index=(" << result.mi << ", " << result.uei << ", " << result.umui << ")"
		  << std::endl;
	iresult++;
      }
      std::cout << "Number of sorts required: " << num_sorts << std::endl;
      auto stop_gridscan = std::chrono::high_resolution_clock::now();
      auto duration_gridscan = std::chrono::duration_cast<std::chrono::milliseconds>(stop_gridscan - start_gridscan);  
  
      // we use the top grid points to seed a simplex search
      //minimizer.algoName = "Migrad";
      minimizer.algoName = "Simplex";

      // over thinking this
      // double seed_par[3] = {       
      // if ( reduce_dim>=0 ) {
      // 	if ( reduce_method=="marginalize" ) {
      // 	  minimizer.fixOnePar( reduce_dim, xxx );
      // 	}
      // 	else if ( reduce_method=="profile" ) {
      // 	  minimizer.onlyFitOnePar( reduce_dim, xxx );
      // 	}
      // }

      // trying a variety of starts
      std::vector< simplex_result_t > simplex_v;
  
      double min_minimizer =100000000;
      std::vector<double> beststart(4,1.0e3);
      
      auto start_simplex = std::chrono::high_resolution_clock::now();
      
      if(DO_SIMPLEX_SEARCH){
	
	std::vector<double> temp;
	
	std::cout<<"starting data fit loop"<<std::endl;
	
	for (int iiresult=0; iiresult<NUM_TOP_FOR_SIMPLEX; iiresult++) {
	  if ( iiresult>=result_list_v.size() )
	    break;
	  
	  auto const& result = result_list_v.at(iiresult);      

	  int uei_base  = result.uei;
	  int umui_base = result.umui;
	  
	  float mnu_base = sqrt( griddef.dm2_bin_v[ result.mi ] );
	  float ue_base  = griddef.Ue4_bin_v[uei_base];
	  float um_base  = griddef.Umu4_bin_v[umui_base];
	  
	  temp = minimizer.doFit( data_v, mnu_base, ue_base , um_base );
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
      }//end of do simplex
  
      std::sort( simplex_v.begin(), simplex_v.end() );
  
      std::cout << "SIMPLEX RESULTS" << std::endl;
      int isimplex = 0;
      for (auto& simplex : simplex_v ) {
	std::cout << "[" << isimplex << "] chi2=" << simplex.llh
		  << " mu=" << simplex.mi << " Ue4=" << simplex.uei << " Umu4=" << simplex.umui
		  << std::endl;
      }
      auto stop_simplex = std::chrono::high_resolution_clock::now();
      auto duration_simplex = std::chrono::duration_cast<std::chrono::milliseconds>(stop_simplex - start_simplex);  
  
      // FINAL MIGRAD FIT
      auto start_migrad = std::chrono::high_resolution_clock::now();  
      minimizer.algoName = "MIGRAD";
      std::vector<double> final_result = minimizer.doFit( data_v, simplex_v.front().mi, simplex_v.front().uei, simplex_v.front().umui );
      std::cout << "FINAL MIGRAD RESULT" << std::endl;
      std::cout << "[0] chi2=" << final_result[0]
		<< " mu=" << final_result[1] << " Ue4=" << final_result[2] << " Umu4=" << final_result[3]
		<< std::endl;

      auto stop_migrad = std::chrono::high_resolution_clock::now();
      auto duration_migrad = std::chrono::duration_cast<std::chrono::milliseconds>(stop_migrad - start_migrad);  
      
      // SAVE BEST FIT
      if ( simplex_v.front().llh < final_result[0] ) {
	// save simplex output
	beststart[0] = simplex_v.front().llh;
	beststart[1] = simplex_v.front().mi*simplex_v.front().mi;
	beststart[2] = simplex_v.front().uei;
	beststart[3] = simplex_v.front().umui;
      }
      else  {
	beststart = final_result;
      }
      
      // Recalc final chi2
      double final_bestfit_par[3] = { log(beststart[1]), beststart[2], beststart[3] };
      double final_chi2 = sbn::SBNllminimizer::negative_likelihood_ratio( final_bestfit_par );
      
      if ( fabs(beststart[0]-final_chi2)>1e-5 ) {
	std::cout << "FINAL CHI2 and BEST CHI2 DISAGREES!" << std::endl;
	std::cout << "  beststart chi2 = " << beststart[0] << std::endl;
	std::cout << "  final chi2 = " << final_chi2 << std::endl;
      }
  
      // chifile<<beststart[0]<<" "<<beststart[1]<<" "<<beststart[2]<<" "<<beststart[3] << " "
      // 	     << sbn::SBNllminimizer::llh_components_v[1] << " " << sbn::SBNllminimizer::llh_components_v[2]
      // 	     << std::endl;

      beststart.push_back( sbn::SBNllminimizer::llh_components_v[1] );
      beststart.push_back( sbn::SBNllminimizer::llh_components_v[2] );
  
      // float ue4_val = 0;
      // float umu4_val = .21;
      // float dm2_val = sqrt(3.37);
      // finally get some stuff to print...
      float e_app_best = 4*pow(beststart[2],2)*pow(beststart[3],2);
      float e_dis_best = 4*pow(beststart[2],2)*(1-pow(beststart[2],2));
      float m_dis_best = 4*pow(beststart[3],2)*(1-pow(beststart[3],2));
      std::cout<<e_app_best<<" "<<e_dis_best<<" "<<m_dis_best<<std::endl;

      std::cout << "DURATION GRIDSCAN: " << duration_gridscan.count() << std::endl;
      std::cout << "DURATION SIMPLEX: "  << duration_simplex.count() << std::endl;
      std::cout << "DURATION MIGRAD: "   << duration_migrad.count() << std::endl;
      
      return beststart;
      
    };//end of analyze procedure
    
  };//end of class definition



}//end of namespace

#endif
