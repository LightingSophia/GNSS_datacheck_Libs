
#ifndef GXTRGRC_H
#define GXTRGRC_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
 
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements extration class
  Version: $ Rev: $

  2013-01-10 /JD: created

-*/

#include "../gio/gxtrqc.h"

using namespace std;
using namespace pugi;

namespace gnut {
  
class t_gxtrgrc : public t_gxtrqc
{
 public:
   t_gxtrgrc(t_gsetbase* set, const string& pgm, const t_gtime& dt);
   t_gxtrgrc();   
  
   virtual ~t_gxtrgrc();

 protected:
   virtual void _get_settings();
   
   virtual void _grckpi(ostringstream& os);      
   virtual void _grchis(ostringstream& os);

   t_gtriple _pos_ref;
   
   double _max_vpe;
   double _max_hpe;
   double _max_dop;
         
   double _hist_min;
   double _hist_max;
   double _hist_int;
};

} // namespace

#endif
