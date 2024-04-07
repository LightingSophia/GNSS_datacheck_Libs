
#ifndef GSETPROC_H
#define GSETPROC_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements processing setting
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#define XMLKEY_PROC "process"

#include <string>
#include <iostream>

#include "../gio/glog.h"
#include "../gutils/gtypeconv.h"
#include "../gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut {

enum CONSTRPAR { EST, FIX };      
   
enum GRDMPFUNC { DEF_GRDMPFUNC, TILTING, CHEN_HERRING, BAR_SEVER };
enum ZTDMPFUNC { DEF_ZTDMPFUNC, COSZ, GMF };
enum IONMPFUNC { DEF_IONMPFUNC, ICOSZ, QFAC, NONE };
enum OBSWEIGHT { DEF_OBSWEIGHT, EQUAL, SINEL, SINEL2, SINEL4, CODPHA, MLTPTH };
enum TROPMODEL { DEF_TROPMODEL, SAASTAMOINEN, DAVIS, HOPFIELD, MOPS, GPTW, GPT2W, GAL27, GALTROPO27, EXTERN };
enum RESIDTYPE { DEF_RESIDTYPE, RES_ORIG, RES_NORM, RES_ALL };
enum OBSCOMBIN { DEF_OBSCOMBIN, IONO_FREE, RAW_SINGLE, RAW_DOUBLE, RAW_ALL };
enum PHASEBIAS { DEF_POSTPRD,  RTCM_CORR, ICLK_CORR};
enum ATTITUDES { DEF_YAWMODEL, YAW_NOMI, YAW_RTCM };
enum CBIASCHAR { DEF_CBIASCHAR, ORIG, CHAR2, CHAR3};
   
class t_gsetproc : public virtual t_gsetbase
{
 public:
   t_gsetproc();
  ~t_gsetproc();
   
   void check();                                 // settings check
   void help();                                  // settings help

   bool   tropo();
   bool   iono();
   bool   tropo_slant();
   bool   tropo_grad();
   bool   phase();
   bool   pos_kin();
   bool   use_eclipsed();
   bool   clean_iono();
   bool   auto_band();
   double sig_init_ztd();
   double sig_init_vion();
   double sig_init_grd();   
   double sig_init_crd();
   double sig_init_amb();
   double sig_init_glo();
   double sig_init_gal();
   double sig_init_bds();
   double sig_init_qzs();   
   double minimum_elev();
   double max_res_norm();
   
   TROPMODEL tropo_model();   
   ZTDMPFUNC tropo_mf();
   IONMPFUNC iono_mf();   
   GRDMPFUNC grad_mf();   
   CONSTRPAR crd_est();
   OBSWEIGHT weighting();
   RESIDTYPE residuals();
   OBSCOMBIN obs_combin();
   ATTITUDES attitudes();
   CBIASCHAR cbiaschar();

   string    trop_corr();
   string    code_bias_src();
   string    phase_bias_src();   
   
   GRDMPFUNC str2grdmpfunc(string mf);
   ZTDMPFUNC str2ztdmpfunc(string mf);
   IONMPFUNC str2ionmpfunc(string mf);   
   OBSWEIGHT str2obsweight(string wg);
   TROPMODEL str2tropmodel(string tm);
   RESIDTYPE str2residtype(string rs);
   OBSCOMBIN str2obscombin(string oc);
   ATTITUDES str2attitudes(string at);
   CBIASCHAR str2cbiaschar(string cb);   

   string grdmpfunc2str(GRDMPFUNC MF);
   string ztdmpfunc2str(ZTDMPFUNC MF);
   string ionmpfunc2str(IONMPFUNC MF);   
   string obsweight2str(OBSWEIGHT WG);
   string tropmodel2str(TROPMODEL TM);
   string residtype2str(RESIDTYPE RS);
   string obscombin2str(OBSCOMBIN OC);
   string cbiaschar2str(CBIASCHAR CB);   
   string attitudes2str(ATTITUDES AT);
   
 protected:
   
   bool   _phase;                                // use phase data [true=1, false=0]
   bool   _tropo;                                // estimate troposphere [true=1, false=0]   
   bool   _iono;                                 // estimate ionosphere [true=1, false=0]   
   bool   _tropo_grad;                           // estimate troposphere gradinet [true=1, false=0]
   bool   _tropo_slant;                          // estimate tropo slant delays
   TROPMODEL _tropo_model;                       // tropospheric model [SAASTAMOINEN, DAVIS, HOPFIELD, ...]
   ZTDMPFUNC _tropo_mf;                          // tropo mapping function [COSZ, NIELL, GMF, VMF1, ... ]
   IONMPFUNC _iono_mf;                           // iono mapping function [COSZ, QFAC, NONE, ... ]   
   GRDMPFUNC _grad_mf;                           // grad mapping function [tilt, CHH]
   OBSWEIGHT _obs_weight;                        // observations weighting
   RESIDTYPE _res_type;                          // types of residuals
   OBSCOMBIN _obs_combin;                        // observation combination 
   ATTITUDES _attitudes;						             // satellite attitude model
   CBIASCHAR _cbiaschar;						             // forcing code bias signals

   double _sig_init_ztd;                         // accuracy of initial zenith path delays [m]
   double _sig_init_vion;                        // accuracy of initial vertical iono path delays [m]
   double _sig_init_grd;                         // accuracy of initial tropo gradients [m]   
   double _sig_init_crd;                         // accuracy of initial coordinates [m]
   double _sig_init_amb;                         // accuracy of initial ambiguity [m]
   double _sig_init_glo;                         // accuracy of initial GLONASS system time difference
   double _sig_init_gal;                         // accuracy of initial Galileo system time difference
   double _sig_init_bds;                         // accuracy of initial BeiDou system time difference
   double _sig_init_qzs;                         // accuracy of initial QZSS system time difference
   double _minimum_elev;                         // elevation angle cut-off [degree]
   double _max_res_norm;                         // normalized residuals threshold
   string _crd_est;                              // FIX or estimate CRD
   bool   _pos_kin;                              // static/kinematic receiver (true == kinematic)
   bool   _use_eclipsed;                         // used eclipsed satellites
   bool   _auto_band;                            // automatic band selection (mainly for Anubis purpose)
   bool   _clean_iono;                           // estimated rec DCB for cleaning ionosphere estimates
   string _code_bias_src;                        // source for code  bias
   string _phase_bias_src;                       // source for phase bias
   
 private:
};

} // namespace

#endif
