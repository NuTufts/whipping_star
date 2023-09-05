#ifndef __GRID_DEFINITION_H__
#define __GRID_DEFINITION_H__

#include <vector>
#include "TMath.h"

class GridDefinition {

 public:

  GridDefinition() {};
  virtual ~GridDefinition() {};

  std::vector<float> dm2_bin_v;
  std::vector<float> Ue4_bin_v;
  std::vector<float> Umu4_bin_v;
  std::vector< std::vector<float> > oscpar_v;
  std::vector< std::vector<int> >   binindex_v;


  double dm2_lowbound;
  double dm2_hibound;
  double ue4_lowbound;
  double ue4_hibound;
  double umu4_lowbound;
  double umu4_hibound;

  int dm2_grdpts;
  int ue4_grdpts;
  int umu4_grdpts;

  void clear() {
    dm2_bin_v.clear();
    Ue4_bin_v.clear();
    Umu4_bin_v.clear();
  };

  int binindex_to_kindex( int mi, int uei, int umui ) {
    return mi*ue4_grdpts*umu4_grdpts + uei*umu4_grdpts + umui;
  };
  
  void define_maya_grid() {

    // Grid definition
    
    dm2_grdpts  = 100;
    ue4_grdpts  = 50;
    umu4_grdpts = 50;

    ue4_lowbound  = 0.01;
    ue4_hibound   = 0.5;
    //ue4_hibound   = 0.70710678118;
    umu4_lowbound = 0.01;
    umu4_hibound  = 0.5;
    //umu4_hibound  = 0.70710678118;

    oscpar_v.reserve( dm2_grdpts*ue4_grdpts*umu4_grdpts );
    binindex_v.reserve( dm2_grdpts*ue4_grdpts*umu4_grdpts );
    
    // Bin Edge definitions
    
    const float dm2_maya[101] = {0.0101158, 0.0110917, 0.0121619, 0.0133352, 0.0146218, 0.0160324, 0.0175792, 0.0192752, 0.0211349, 0.0231739,
				 0.0254097, 0.0278612, 0.0305492, 0.0334965, 0.0367282, 0.0402716, 0.044157, 0.0484172, 0.0530884, 0.0582102,
				 0.0638262, 0.0699841, 0.076736, 0.0841393, 0.0922569, 0.101158, 0.110917, 0.121618, 0.133352, 0.146217,
				 0.160324, 0.175792, 0.192752, 0.211348, 0.231739, 0.254096, 0.278611, 0.305491, 0.334964, 0.367281, 0.402716,
				 0.441569, 0.484171, 0.530882, 0.582101, 0.638261, 0.699839, 0.767358, 0.841392, 0.922567, 1.01158, 1.10917,
				 1.21618, 1.33352, 1.46217, 1.60324, 1.75791, 1.92752, 2.11348, 2.31738, 2.54096, 2.78611, 3.0549, 3.34964,
				 3.6728, 4.02715, 4.41568, 4.8417, 5.30881, 5.821, 6.3826, 6.99838, 7.67357, 8.4139, 9.22565, 10.1157, 11.0917,
				 12.1618, 13.3351, 14.6217, 16.0323, 17.5791, 19.2751, 21.1347, 23.1738, 25.4095, 27.861, 30.549, 33.4963, 36.7279,
				 40.2714, 44.1567, 48.4168, 53.088, 58.2098, 63.8258, 69.9836, 76.7355, 84.1388, 92.2563, 112.202};

    dm2_lowbound = dm2_maya[0];
    dm2_hibound  = dm2_maya[99];
    
    dm2_bin_v.resize(dm2_grdpts+1,0);
    for (int i=0; i<=dm2_grdpts; i++) {
      dm2_bin_v[i] = dm2_maya[i];
    }

    Ue4_bin_v.resize(ue4_grdpts+1,0.0);
    for (int uei_base=0; uei_base<ue4_grdpts; uei_base++) {
      //float ue_base = pow(10.,(uei_base/float(ue4_grdpts-1)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
      float ue_base = pow(10.,(uei_base/float(ue4_grdpts)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
      Ue4_bin_v[uei_base] = ue_base;
    }

    Umu4_bin_v.resize(umu4_grdpts+1,0.0);
    for (int umui_base=0; umui_base<umu4_grdpts; umui_base++) {
      //float umu_base = pow(10.,(umui_base/float(umu4_grdpts-1)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
      float umu_base = pow(10.,(umui_base/float(umu4_grdpts)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
      Umu4_bin_v[umui_base] = umu_base;
    }

    // define list of oscpars
    int kindex = 0;
    for(int mi_base = 0; mi_base < dm2_grdpts; mi_base += 1 ){
      for(int uei_base = 0; uei_base < ue4_grdpts; uei_base += 1 ){
	for(int umui_base = 0; umui_base < umu4_grdpts; umui_base += 1 ) {
	  std::vector<float> gridpt_oscpar = { dm2_bin_v[mi_base], Ue4_bin_v[uei_base], Umu4_bin_v[umui_base] };
	  oscpar_v.push_back( gridpt_oscpar );
	  std::vector<int>   gridpt_index  = { mi_base, uei_base, umui_base };
	  binindex_v.push_back( gridpt_index );
	  kindex++;
	}
      }
    }
    return;
  };


  void define_proper_bounds_course() {

    // Grid definition    
    dm2_grdpts  = 25;
    ue4_grdpts  = 10;
    umu4_grdpts = 10;

    dm2_lowbound  = 0.01;
    dm2_hibound   = 100.0;
    ue4_lowbound  = 0.01;
    ue4_hibound   = 0.70710678118;
    umu4_lowbound = 0.01;
    umu4_hibound  = 0.70710678118;    

    oscpar_v.reserve( dm2_grdpts*ue4_grdpts*umu4_grdpts );
    binindex_v.reserve( dm2_grdpts*ue4_grdpts*umu4_grdpts );
    
    dm2_bin_v.resize(dm2_grdpts,0);
    for (int i=0; i<dm2_grdpts; i++) {
      float mnu_base = pow(10.0, (i/float(dm2_grdpts-1))*TMath::Log10(dm2_hibound/dm2_lowbound) + TMath::Log10(dm2_lowbound) );
      dm2_bin_v[i] = mnu_base;
    }
    
    Ue4_bin_v.resize(ue4_grdpts,0.0);
    for (int uei_base=0; uei_base<ue4_grdpts; uei_base++) {
      float ue_base = pow(10.,(uei_base/float(ue4_grdpts-1)*TMath::Log10(ue4_hibound/ue4_lowbound) + TMath::Log10(ue4_lowbound)));
      Ue4_bin_v[uei_base] = ue_base;
    }
    
    Umu4_bin_v.resize(umu4_grdpts,0.0);
    for (int umui_base=0; umui_base<umu4_grdpts; umui_base++) {
      float umu_base = pow(10.,(umui_base/float(umu4_grdpts-1)*TMath::Log10(umu4_hibound/umu4_lowbound) + TMath::Log10(umu4_lowbound)));
      Umu4_bin_v[umui_base] = umu_base;
    }
    
    // define list of oscpars
    int kindex = 0;
    for(int mi_base = 0; mi_base < dm2_grdpts; mi_base += 1 ){
      for(int uei_base = 0; uei_base < ue4_grdpts; uei_base += 1 ){
	for(int umui_base = 0; umui_base < umu4_grdpts; umui_base += 1 ) {
	  std::vector<float> gridpt_oscpar = { dm2_bin_v[mi_base], Ue4_bin_v[uei_base], Umu4_bin_v[umui_base] };
	  oscpar_v.push_back( gridpt_oscpar );
	  std::vector<int>   gridpt_index  = { mi_base, uei_base, umui_base };
	  binindex_v.push_back( gridpt_index );
	  kindex++;
	}
      }
    }
    return;
  };
  
};

#endif
