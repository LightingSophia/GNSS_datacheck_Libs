
#ifndef GSTAT_H
#define GSTAT_H 
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: statistical function (1D)
  Version: $ Rev: $

  2014-04-18 /PV: created
  2018-09-28 /JD: revised

-*/

#include "../gdata/gdata.h"
#include "../gutils/gpair.h"

#include <map>
#include <set>
#include <vector>

#define CONF_INT 3    // confident interval factor

using namespace std;

namespace gnut {

class t_gstat : public t_gdata 
{
 public:
   
   typedef map<t_gpair, int>  t_hist;

   t_gstat(double cint = CONF_INT);                           // just construct
   t_gstat(vector<double>& data, double cint = CONF_INT);     // construct + calculate statistics
   virtual ~t_gstat();

   void add_data(vector<double>& data);                       // clear + calculate statistics

   virtual bool valid(){ return _valid; }                     // get status

   virtual int calc_stat(double sig = 0.0);                   // calculate statistics

   virtual int calc_quartiles( double& low, double& upp );
   virtual int calc_iqrlimits( double& low, double& upp );

// virtual double calc_mad(){    _calc_mad();    return _mad;    }

   virtual double get_min(){                     return _min;    }
   virtual double get_max(){                     return _max;    }
   virtual double get_var(){                     return _var;    }
   virtual double get_rms(){                     return _rms;    }
   virtual double get_sdev(){                    return _sdev;   }
   virtual double get_mean(){                    return _mean;   }
   virtual double get_median(){                  return _medi;   }
// virtual double get_mad(){                     return _mad;    }
   
   virtual int    get_size(){                    return _data.size(); }
   virtual int    get_outl(){                    return _outl.size(); }

   virtual vector<double>  get_outl_data(){      return _outl; }   

   t_hist histogram(vector<double>& data, set<double>& bound);

   // temporary solution for weighted mean
   double calc_wmean(vector<pair<double,double>> data);

 protected:
   
   void         _add_data(vector<double>& data);

   virtual void _clear();

   // calculate statistics with outlier detection and set validity status
   virtual int  _statistics(double sig = 0.0);
   
   // internal purposes only (no validity status changed)
   virtual int  _calc_minmax();
   virtual int  _calc_median();
   virtual int  _calc_mean();
   virtual int  _calc_sdev();
   virtual int  _calc_rms();
// virtual int  _calc_mad();

   virtual int  _chk_resid(double sig = 0.0); // check residuals (optional use of external sigma)

// virtual double _p(double v);

   bool _valid;

   vector<double> _data;
   vector<double> _outl;

   double _cint;
   double _sdev;
   double _var;
   double _rms;
   double _mean;
   double _medi;
// double _mad;
   double _min;
   double _max;
};
   
} // namespace

#endif
