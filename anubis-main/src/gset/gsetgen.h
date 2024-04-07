
#ifndef GSETGEN_H
#define GSETGEN_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements common general settings
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#define XMLKEY_GEN   "gen"

#define DEF_RECEIVER "   " // all !
#define DEF_SAMPLING 30    // 30s !

#include <set>
#include <string>
#include <iostream>

#include "../gio/glog.h"
#include "../gutils/gsys.h"
#include "../gutils/gtime.h"
#include "../gutils/gtypeconv.h"
#include "gsetbase.h"

using namespace std;

namespace gnut {

class t_gsetgen : public virtual t_gsetbase
{
 public:
   t_gsetgen( bool gnss = true );                 // distinguish GNSS/nonGNSS app
  ~t_gsetgen();

   void check();                                  // settings check
   void help();                                   // settings help
   
   t_gtime beg();
   t_gtime end();
   double sampling();
   double sampling_default(){ return DEF_SAMPLING; }       // default sampling
   int    sampling_scalefc(){ return (int)pow(10,_dec); }  // decimals scale
   int    sampling_decimal(){ return _dec; }               // decimals for sampling interval (for high-rate)

   virtual set<string> sys();
   virtual set<string> rec();

 protected:

   bool   _gnss;
   string _sys;
   int    _dec;
      
 private:
};

} // namespace

#endif
