
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
 
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.
 
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <cmath>

#include "../gio/gxtrgrc.h"
#include "../gset/gsetqc.h"
#include "../gutils/gstat.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructors
// ---------
t_gxtrgrc::t_gxtrgrc()
{
}
   
t_gxtrgrc::t_gxtrgrc(t_gsetbase* set,
             const string& pgm, const t_gtime& dt)
 : t_gxtrqc(set, pgm, dt)
{
  gtrace("t_gxtrgrc::t_gxtrgrc");

  _get_settings();      
  _qcdata = make_shared<t_gqcdata>();

  _hist_min = -4.0;
  _hist_max =  4.0;
  _hist_int =  0.5;

}

// Destructor
// ---------
t_gxtrgrc::~t_gxtrgrc()
{
  gtrace("t_gxtrgrc::desctructor");
}


// ===================
// PROTECTED FUNCTIONS
// ===================
   
// Get settings
// -----------------------------
void t_gxtrgrc::_get_settings()
{
  gtrace("t_gxtrgrc::_get_settings");

  t_gxtrqc::_get_settings();
   
  _gkpi             = dynamic_cast<t_gsetqc*>(_set)->gkpi();
  _hist             = 0;
  
  _max_vpe          = dynamic_cast<t_gsetqc*>(_set)->max_vpe();
  _max_hpe          = dynamic_cast<t_gsetqc*>(_set)->max_hpe();
  _max_dop          = dynamic_cast<t_gsetqc*>(_set)->max_dop();

  if(  _gkpi && !_band ) _band = -1; // do always band testing
  if(  _gkpi && !_calc ) _calc = -1; // do always position calculation
  if(  _gkpi && !_summ ) _summ = -1; // maybe separate in the future (counts of epoch)

}

   
// Section: KPI pars
// -----------------
void t_gxtrgrc::_grckpi(ostringstream& os)
{
  gtrace("t_gxtrgrc::_grckpi");

  if( _gkpi > 0 ){ _section(os, "Key-parameter indicators", _gkpi);
  }else if( !_gkpi ) return;

  set<GSYS> set_sys = _gnav->systems();

  string tt(_ref.str_ymdhms() + " ");
  
  ostringstream os_v2;  
  
  // loop GSYS
  for(auto itSYS = set_sys.begin(); itSYS != set_sys.end(); ++itSYS ){
    string gs = t_gsys::gsys2str(*itSYS);
//  char   ch = t_gsys::gsys2char(*itSYS);

    // SETTINGS INFORMATION
    _pos_ref = _qcdata->qc_est.xyz_crd[*itSYS];
    shared_ptr<t_grec> setrec = _rec->grec(_site);
    if( setrec == 0 ){
       if( _log ) _log->comment(0,"gxtrgrc","Error: KPI no user receiver settings "+_site);
       else               cerr << "gxtrgrc - Error: KPI no user receiver settings "+_site << endl;
       setrec = make_shared<t_grec>();      
    }
    else{ _pos_ref = setrec->crd(_obs_beg); }
    
    t_gtriple _blh_ref;
    xyz2ell(_pos_ref, _blh_ref, false);

    string key = gs + "KPI";

    if( _gkpi >= 2 ){
      _setkey(os_v2, key, '#', tt);
      os_v2 << W(9)  << "dN[m]" <<  W(9) << "dE[m]"  << W(9) << "dU[m]"
            << W(9)  << "PDOP"  <<  W(9) << "FLAGS"  << endl;
    }

    for(auto itEPO  = _qcdata->qc_est.pos_line[*itSYS].begin();
             itEPO != _qcdata->qc_est.pos_line[*itSYS].end();
           ++itEPO )
    {
      t_gtriple dNEU;
      t_gtriple dXYZ( itEPO->second.xyz - _pos_ref );

      xyz2neu( itEPO->second.blh, dXYZ, dNEU);
      _qcdata->qc_kpi.nepo[*itSYS]++;

      double pd = itEPO->second.gdop;
      double dh = sqrt(dNEU[0]*dNEU[0] + dNEU[1]*dNEU[1]);
      double dv = abs(dNEU[2]);
      
      string flag  = ( pd > _max_dop ) ? "D" : "_"; // individual flags
             flag += ( dh > _max_hpe ) ? "H" : "_"; // individual flags 
             flag += ( dv > _max_vpe ) ? "V" : "_"; // individual flags

      // conditional summing for statistics
      if( pd  > _max_dop ){   _qcdata->qc_kpi.PD_XX[*itSYS]++;
      }else{                  _qcdata->qc_kpi.PD_OK[*itSYS]++;

  	if(  dh > _max_hpe ){ _qcdata->qc_kpi.dH_XX[*itSYS]++; }
        else{                 _qcdata->qc_kpi.dH_OK[*itSYS]++; }

	if(  dv > _max_vpe ){ _qcdata->qc_kpi.dV_XX[*itSYS]++; }
        else{                 _qcdata->qc_kpi.dV_OK[*itSYS]++; }
      }

        
      if( _gkpi >= 2 ){
        _setkey(os_v2, key, ' ', itEPO->first.str_ymdhms());
        os_v2 << " " << fixed
              << setprecision(3) << W(9) << dNEU[0]
                                 << W(9) << dNEU[1]
                                 << W(9) << dNEU[2]
              << setprecision(2) << W(9) << _qcdata->qc_est.pos_line[*itSYS][itEPO->first].gdop
                                 << W(9) << flag
              << endl;
      }
    }

    ostringstream ostr;
    _setkey(os, "GNSKPI", '#', tt);
    os << W(9) << "CountOK" << W(9) << "CountXX" << W(9) << "%_Expect" << W(9) << "%_Exist" << W(9) << "%_Posit" << W(9) << "%_DopOK" << endl;
 
    int expEpo = _nEpoExpect[*itSYS];
    int havEpo = _nEpoExists[*itSYS];
    int posEpo = _qcdata->qc_kpi.PD_OK[*itSYS] + _qcdata->qc_kpi.PD_XX[*itSYS]; // # of estimable positions
    int dopEpo = _qcdata->qc_kpi.PD_OK[*itSYS];                                 // # of position with DOP<crit

    if( _nEpoExpect[*itSYS] == 0 || _nEpoExpect[*itSYS] )
    {
      double Vsig =          _qcdata->qc_est.blh_rep[*itSYS][2];
      double Hsig = sqrt(pow(_qcdata->qc_est.blh_rep[*itSYS][0],2.0)
                       + pow(_qcdata->qc_est.blh_rep[*itSYS][1],2.0) );

      _setkey(os, gs+"_EP", ' ', tt);
      os << fixed << setprecision(0)
         << W(9) << _qcdata->qc_est.epo_used[*itSYS]
         << W(9) << _qcdata->qc_est.epo_excl[*itSYS]
         << W(9) << dbl2str( (expEpo) ? _qcdata->qc_est.epo_used[*itSYS]*100.0 / expEpo : 0.0 )
         << W(9) << dbl2str( (havEpo) ? _qcdata->qc_est.epo_used[*itSYS]*100.0 / havEpo : 0.0 )
         << W(9) << "-"
         << W(9) << "-" << endl;

      _setkey(os, gs+"_DH", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.dH_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.dH_XX[*itSYS]
         << W(9) << dbl2str( (expEpo) ? _qcdata->qc_kpi.dH_OK[*itSYS]*100.0 / expEpo : 0.0 )
         << W(9) << dbl2str( (havEpo) ? _qcdata->qc_kpi.dH_OK[*itSYS]*100.0 / havEpo : 0.0 )
         << W(9) << dbl2str( (posEpo) ? _qcdata->qc_kpi.dH_OK[*itSYS]*100.0 / posEpo : 0.0 )
         << W(9) << dbl2str( (dopEpo) ? _qcdata->qc_kpi.dH_OK[*itSYS]*100.0 / dopEpo : 0.0 )
         << W(9) << dbl2str(Hsig*1.0)
         << endl;

      _setkey(os, gs+"_DV", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.dV_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.dV_XX[*itSYS]
         << W(9) << dbl2str( (expEpo) ? _qcdata->qc_kpi.dV_OK[*itSYS]*100.0 / expEpo : 0.0 )
         << W(9) << dbl2str( (havEpo) ? _qcdata->qc_kpi.dV_OK[*itSYS]*100.0 / havEpo : 0.0 )
         << W(9) << dbl2str( (posEpo) ? _qcdata->qc_kpi.dV_OK[*itSYS]*100.0 / posEpo : 0.0 )
         << W(9) << dbl2str( (dopEpo) ? _qcdata->qc_kpi.dV_OK[*itSYS]*100.0 / dopEpo : 0.0 )
         << W(9) << dbl2str(Vsig*1.0)
         << endl;
/*
      _setkey(os, gs+"_HD", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.HD_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.HD_XX[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_kpi.HD_OK[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.HD_OK[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.HD_OK[*itSYS]*100.0 /  posEpo) << endl;
//         << W(9) << dbl2str(_qcdata->qc_kpi.HD_OK[*itSYS]*100.0 /  dopEpo) << endl;
*/
      _setkey(os, gs+"_PD", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.PD_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.PD_XX[*itSYS]
         << W(9) << dbl2str( (expEpo) ? _qcdata->qc_kpi.PD_OK[*itSYS]*100.0 / expEpo : 0.0 )
         << W(9) << dbl2str( (havEpo) ? _qcdata->qc_kpi.PD_OK[*itSYS]*100.0 / havEpo : 0.0 )
         << W(9) << dbl2str( (posEpo) ? _qcdata->qc_kpi.PD_OK[*itSYS]*100.0 / posEpo : 0.0 )
         << W(9) << dbl2str( (dopEpo) ? _qcdata->qc_kpi.PD_OK[*itSYS]*100.0 / dopEpo : 0.0 )
	 << endl;
/*
      _setkey(os, gs+"_VD", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.VD_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.VD_XX[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_kpi.VD_OK[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.VD_OK[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.VD_OK[*itSYS]*100.0 /  posEpo)
         << W(9) << dbl2str(_qcdata->qc_kpi.VD_OK[*itSYS]*100.0 /  dopEpo) << endl;

      _setkey(os, gs+"_GD", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.GD_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.GD_XX[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_kpi.GD_OK[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.GD_OK[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.GD_OK[*itSYS]*100.0 /  posEpo)
         << W(9) << dbl2str(_qcdata->qc_kpi.GD_OK[*itSYS]*100.0 /  dopEpo) << endl;
*/
      _setkey(os, gs+"_DF", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.DF_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.DF_XX[*itSYS]
         << W(9) << dbl2str( (expEpo) ? _qcdata->qc_kpi.DF_OK[*itSYS]*100.0 / expEpo : 0.0 )
	 << W(9) << dbl2str( (havEpo) ? _qcdata->qc_kpi.DF_OK[*itSYS]*100.0 / havEpo : 0.0 )
         << W(9) << "-"
         << W(9) << "-" << endl;

      os << endl;
    }
    else{
      if( _log ) _log->comment(0,"gxtrgrc","Warning: KPI ["+gs+"] skipped, no epoch available"); 
      os                              << "# Warning: KPI ["+gs+"] skipped, no epoch available\n";
    }
  }
  
  if(_gkpi >= 2) os << os_v2.str();

}

// Section: Histogram of position residuals
// -----------------
void t_gxtrgrc::_grchis(ostringstream& os)
{
  gtrace("t_gxtrgrc::_histogram");

  if( _hist > 0 ){ _section(os, "Histogram of position residuals", _hist);
  }else if( !_hist ) return;

  vector<double> resH;  
  set<double>    bound;

  // generation boundaries
  for(double b = _hist_min; b <= _hist_max; b += _hist_int){
    bound.insert(b);
  }
  
  set<GSYS> set_sys = _gnav->systems();

  // loop GSYS
  for(auto itSYS = set_sys.begin(); itSYS != set_sys.end(); ++itSYS ){
    string gs = t_gsys::gsys2str(*itSYS);
    
    for(auto itEPO  = _qcdata->qc_est.pos_line[*itSYS].begin();
        itEPO != _qcdata->qc_est.pos_line[*itSYS].end();
        ++itEPO )
    {
      t_gtriple dNEU;
      t_gtriple dXYZ( itEPO->second.xyz - _pos_ref );
      
      xyz2neu( itEPO->second.blh, dXYZ, dNEU);
      resH.push_back(dNEU[2]);
    }
  }

  t_gstat stat;
  map<t_gpair, int> frequency = stat.histogram(resH, bound);

//#ifdef DEBUG
  cout << "Histogram" << endl << "-----------------------" << endl;
  size_t i = 1;
  for(auto it = frequency.begin(); it != frequency.end(); it++){
    ostringstream lstr, pstr;
    if(it != frequency.begin()) lstr << fixed << setprecision(3) << setw(10) << right << it->first.crd(0) << setw(3) << "< ";
    else lstr << setw(10) << "" << right << setw(3) << "  ";
    
    if(i != frequency.size()) pstr << "< " << fixed << setprecision(3) << setw(10) << left << it->first.crd(1);
    else pstr << "  " << setw(10) << "";

    cout << lstr.str() <<  "x " << pstr.str() << setw(4) << "|" << setprecision(0) << setw(4) << it->second << endl;
    i++;
  }
//#endif
  
}
  
} // namespace
