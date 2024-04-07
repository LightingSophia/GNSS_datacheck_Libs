
#ifndef GQCDATA_H
#define GQCDATA_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: container for QC data
  Version: $Rev:$

  2018-05-20 /JD: created

-*/

#include <map>
#include <vector>

#include "../gdata/gdata.h"
#include "../gdata/gsatdata.h"
#include "../gutils/gtime.h"
#include "../gutils/gnss.h"

using namespace std;

namespace gnut {

struct t_qc_est {
  struct t_pos_line {
    t_gtriple xyz,blh;
    double gdop,pdop,hdop,vdop,clk;
    int nsat, xsat;
  };
  double                sampling;
  t_gtime               epo_beg, epo_end;
  map<GSYS,t_gtriple>   xyz_crd, xyz_rep, blh_crd, blh_rep;
  map<GSYS, int>        epo_used, epo_excl;

  map<GSYS,map<t_gtime,t_pos_line>> pos_line;
};

struct t_qc_kpi {
  map<GSYS,int>         nepo;
  map<GSYS,double>      dH_OK, dV_OK, HD_OK, VD_OK, PD_OK, GD_OK, DF_OK;
  map<GSYS,double>      dH_XX, dV_XX, HD_XX, VD_XX, PD_XX, GD_XX, DF_XX;
};

// time-dependent data
struct t_qc_tim {
  struct t_qctimv {
    double ele, azi;
    int    lbn, cbn;
    int    health;
    map<GOBS, double>           snr;
    map<GOBS, pair<double,int>> mpt;
  };

  map<GSYS, map<t_gtime, map<string, t_qctimv>>> dat;
};

struct t_qc_bnd {
  map<GSYS, map<t_gtime, map<string, pair<int,int>>>> data;
};


class t_gqcdata; typedef shared_ptr<t_gqcdata> t_spt_qcdata;
class t_gqcdata : public t_gdata 
{
 public:
  t_gqcdata();
 ~t_gqcdata();

  t_gtime    qc_epo;
  t_qc_est   qc_est;
  t_qc_tim   qc_tim;
  t_qc_bnd   qc_bnd;
  t_qc_kpi   qc_kpi;

};

} // namespace

#endif
