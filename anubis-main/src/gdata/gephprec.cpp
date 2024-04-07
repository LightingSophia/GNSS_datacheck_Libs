
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "../gdata/gephprec.h"
#include "../gutils/gconst.h"
#include "../gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gephprec::t_gephprec() 
 : t_geph(),
   _degree(0),
   _poly_x(),
   _poly_y(),
   _poly_z(),
   _poly_c() 
{
  id_type(  t_gdata::EPHPREC );
  id_group( t_gdata::GRP_EPHEM );
}

// destructor
// ----------
t_gephprec::~t_gephprec()
{}

// add data
// ---------
int t_gephprec::add( string sat, vector<t_gtime> t,
		     const vector<double>& x, const vector<double>& y,
	             const vector<double>& z, const vector<double>& c )
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  // check dimensions
  if( ndata() != t.size() || ndata() != x.size() ||
      ndata() != y.size() || ndata() != z.size() ||
      ndata() != c.size() || sat.empty() ){ _gmutex.unlock(); return -1; }

  _clear();
  _sat = sat;
  _epoch = t[0] + t[_degree].diff( t[0] )/2; // reference epoch always in the mid of interval!!!
   
  _xcrd = x; //  _xcrd.assign(x[0],x[_degree]);
  _ycrd = y; //  _ycrd.assign(y[0],y[_degree]);
  _zcrd = z; //  _zcrd.assign(z[0],z[_degree]);
  _clkc = c; //  _clkc.assign(c[0],c[_degree]);


  for( unsigned int i = 0; i < ndata(); ++i ){
    _dt.push_back( t[i] - _epoch );
     
#ifdef DEBUG
    cout << _sat
	 << " t = "
         << fixed    
         << setprecision(3)
	 << setw(14) << _epoch.str("%Y-%m-%d %H:%M:%S")
         << setw(14) << _dt[i]
	 << setw(14) << _xcrd[i] 
	 << setw(14) << _ycrd[i] 
	 << setw(14) << _zcrd[i] 
         << setprecision(6)
	 << setw(14) << _clkc[i]*1000000
         << endl;
#endif

  }

  _poly_x.polynomials(_dt,_xcrd);
  _poly_y.polynomials(_dt,_ycrd);
  _poly_z.polynomials(_dt,_zcrd);
  _poly_c.polynomials(_dt,_clkc);

  _gmutex.unlock(); return 0;
}


// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// ----------
int t_gephprec::pos( const t_gtime& t, double xyz[], double var[], double vel[], bool chk_health )
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

            xyz[0] = xyz[1] = xyz[2] = 0.0;
  if( var ) var[0] = var[1] = var[2] = 0.0;
  if( vel ) vel[0] = vel[1] = vel[2] = 0.0;

  double tdiff = t - _epoch;

  if( ! _valid_crd() || sat().empty() ||
      fabs(_epoch-t) > (fabs( _poly_x.xref() - _poly_x.span()/2 ) + 0.25 )  ){
    if( _log ) _log->comment(0, "gephprec:pos", "warning: no ephemeris ["
			         + _sat + t.str("] %Y-%m-%d %H:%M:%S")
			         + " epoch: " + _epoch.str("%Y-%m-%d %H:%M:%S")
			         + " tdiff: " + dbl2str(tdiff)
			         + " xref: " + dbl2str(_poly_x.xref())
			         + " span: " + dbl2str(_poly_x.span())
			        );
    _gmutex.unlock(); return -1;
  }
   
  _poly_x.evaluate(tdiff, 0, xyz[0]);                // earth-fixed coordinates [m]
  _poly_y.evaluate(tdiff, 0, xyz[1]);
  _poly_z.evaluate(tdiff, 0, xyz[2]);

  if( vel ){
    _poly_x.evaluate(tdiff, 1, vel[0]);
    _poly_y.evaluate(tdiff, 1, vel[1]);
    _poly_z.evaluate(tdiff, 1, vel[2]);
  }

  _gmutex.unlock(); return 1;
}
   

// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// NOT CORRECTED FOR 2nd relativistic effect !!!
// ----------
int t_gephprec::clk( const t_gtime& t, double*    clk, double*   var, double*  dclk, bool chk_health )
{    
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

             *clk  = 0.0;
  if( var  ) *var  = 0.0;
  if( dclk ) *dclk = 0.0;

  if( ! _valid_clk() || sat().empty() ){ _gmutex.unlock(); return -1; }  // check if valid
   
  double tdiff = t - _epoch; 

  if (abs(tdiff) > (abs(_poly_c.xref() - _poly_c.span() / 2) + 0.25)){
    if( _log ) _log->comment(0, "gephprec:clk", "warning: no clock [" + _sat + t.str("] %Y-%m-%d %H:%M:%S"));
    _gmutex.unlock(); return -1;
  }
   
  _poly_c.evaluate(tdiff, 0, *clk);

  if( var ){} // not implemented

  if( dclk ){
    _poly_c.evaluate(tdiff, 1, *dclk);
  }
   
  _gmutex.unlock(); return 1;
}

   
// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// ----------
int t_gephprec::pos_int( const t_gtime& t, double xyz[], double var[], double vel[] )
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
            xyz[0] = xyz[1] = xyz[2] = 0.0;
  if( var ) var[0] = var[1] = var[2] = 0.0;
  if( vel ) vel[0] = vel[1] = vel[2] = 0.0;

  if( ! _valid_crd() || sat().empty() ){ _gmutex.unlock(); return -1; } // check if valid

  double tdiff = t - _epoch;

  if( abs(tdiff) > (abs( _poly_x.xref() - _poly_x.span()/2 +0.25)) ){
    if( _log ) _log->comment(0, "gephprec:pos_int", "warning: no ephemeris [" + _sat + t.str("] %Y-%m-%d %H:%M:%S"));
    _gmutex.unlock(); return -1;
  }
   
  for( unsigned int i = 0; i < ndata(); ++i ){
     
#ifdef DEBUG
    cout << _sat << t.str(" %Y-%m-%d %H:%M:%S[%T] ")
	 << " t_int = "
	 << tdiff
         << fixed    
         << setprecision(3)
         << setw(14) << _dt[i]
	 << setw(14) << _xcrd[i] 
	 << setw(14) << _ycrd[i] 
	 << setw(14) << _zcrd[i] 
         << setprecision(6)
	 << setw(14) << _clkc[i]*1000000
         << endl;
#endif
 
  }
   
  t_gpoly poly;
  if( var ){
    poly.interpolate( _dt, _xcrd, tdiff, xyz[0], var[0] );
    poly.interpolate( _dt, _ycrd, tdiff, xyz[1], var[1] );
    poly.interpolate( _dt, _zcrd, tdiff, xyz[2], var[2] );
  }else{	
    double dummy[3];
    poly.interpolate( _dt, _xcrd, tdiff, xyz[0], dummy[0] );
    poly.interpolate( _dt, _ycrd, tdiff, xyz[1], dummy[1] );
    poly.interpolate( _dt, _zcrd, tdiff, xyz[2], dummy[2] );
  }
   
  if( vel ){} // not implemented

  _gmutex.unlock(); return 1;
}

      
// t .. GPS time of transmission,
// i.e. GPS time corrected for transit time (range/speed of light)
// NOT CORRECTED FOR 2nd relativistic effect !!!
// ----------
int t_gephprec::clk_int( const t_gtime& t, double*    clk, double*   var, double*  dclk )
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

             *clk  = 0.0;
  if(  var ) *var  = 0.0;
  if( dclk ) *dclk = 0.0;

  if( ! _valid_clk() || sat().empty() ){ _gmutex.unlock(); return -1; }    // check if valid
   
  double tdiff = t - _epoch;

  if( abs(tdiff) > (abs( _poly_c.xref() - _poly_c.span()/2)+0.25 ) ){
    if( _log ) _log->comment(0, "gephprec:clk_int", "warning: no ephemeris [" + _sat + t.str("] %Y-%m-%d %H:%M:%S"));
    _gmutex.unlock(); return -1;
  }
   
  t_gpoly poly;
   
  if( var ){
    poly.interpolate( _dt, _clkc, tdiff, *clk, *dclk );
  }else{
    double dummy;
    poly.interpolate( _dt, _clkc, tdiff, *clk, dummy );
  }
   
  if( dclk ){} // not yet implemented
     
  _gmutex.unlock(); return 1;
}


// validity
// ----------
bool t_gephprec::valid(const t_gtime& t) const
{
  if( _valid_crd() && ! sat().empty() &&
     fabs(_epoch-t) < fabs( _poly_x.xref() - _poly_x.span()/2 ) + 0.25 )  // add 0.25s!
  {
//    cout << " valid : " << sat() << " " << t.str("%Y-%m-%d %H:%M:%S") << endl;
    return true;
  }

//  cout << " ! not valid : " << sat() << " " << t.str("%Y-%m-%d %H:%M:%S") << endl;
  return false;
}
   


// clean data
// ----------
void t_gephprec::_clear()
{
 _epoch = FIRST_TIME;
 _sat.clear();
   
 _dt.clear();

 _xcrd.clear();
 _ycrd.clear();
 _zcrd.clear();
 _xcrd.clear();
   
 _poly_x.reset();
 _poly_y.reset();
 _poly_z.reset();
 _poly_c.reset();
}


// validity ?
// ----------
bool t_gephprec::_valid_crd() const
{
  if( t_geph::_valid() &&
      _poly_x.valid() &&
      _poly_y.valid() && 
      _poly_z.valid() 
    ) return true;

  return false;
}
   
// validity ?
// ----------
bool t_gephprec::_valid_clk() const
{

  bool undef = false; 
  for(unsigned int i = 0; i < _clkc.size(); i++) {
     if(_clkc[i] >= 999999) {
 	   if( _log ) _log->comment(1, "gephprec", "warning: clk_int undefined ephemeris [" + _sat + _epoch.str("] %Y-%m-%d %H:%M:%S"));
	   undef = true;
     }     
  }
   
  if( t_geph::_valid() &&
      _poly_c.valid() &&
      ! undef
    ) return true;

  return false;
}   

} // namespace
