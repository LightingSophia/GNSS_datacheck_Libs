
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

#include <iostream>
#include <iomanip>

#include "../gutils/gobs.h"

using namespace std;

namespace gnut {

// -------------------------------------------------------------------------------------------
// class T_GATTR
// -------------------------------------------------------------------------------------------
bool t_gattr::valid() const
{
  return ( _gattr != ATTR );
}

// set attr
// ----------
void t_gattr::attr( const GOBSATTR& a )
{
  _gattr = a;
}

// get attr
// -----------
GOBSATTR t_gattr::attr()const
{ 
  return _gattr;
}

// operator
// ----------
bool t_gattr::operator==(const t_gattr& g) const
{
  return ( _gattr == g.attr() );
}


// -------------------------------------------------------------------------------------------
// class T_GBAND
// -------------------------------------------------------------------------------------------
bool t_gband::valid() const
{
  return ( t_gattr::valid() && _gband != BAND );
}

// set band
// ----------
void t_gband::band( const GOBSBAND& b )
{
  _gband = b;
}

// get band
// -----------
GOBSBAND t_gband::band()const
{ 
  return _gband;
}

// set attr
// ----------
void t_gband::gattr( const t_gattr& g )
{
  _gattr = g.attr();
  _gband = BAND;
}

// get attr
// -----------
t_gattr t_gband::gattr()const
{
  t_gattr g(_gattr);
  return g;
}

// operators
// ----------
bool t_gband::operator==(const t_gband& g) const
{
  return ( _gband == g.band() &&
           _gattr == g.attr()
         );
}

// -------------------------------------------------------------------------------------------
// class T_GOBS
// -------------------------------------------------------------------------------------------

// valid ?
// -----------
bool t_gobs::valid() const
{
  return ( t_gband::valid() && _gtype != TYPE );
}

// set type
// ----------
void t_gobs::type( const GOBSTYPE& t )
{
  _gtype = t;
}

// get type
// -----------
GOBSTYPE t_gobs::type()const
{
  return _gtype;
}

// set gband
// ----------
void t_gobs::gband( const t_gband& g )
{
  _gattr = g.attr();
  _gband = g.band();
  _gtype = TYPE;
}

// get gband
// -----------
t_gband t_gobs::gband()const
{
  t_gband g(_gband,_gattr);
  return g;
}

// operator
// ----------
bool t_gobs::operator==(const t_gobs& g) const
{
  return ( _gtype == g.type() &&
           _gband == g.band() &&
           _gattr == g.attr()
   );
}

// set from GOBS
// -----------
int t_gobs::gobs(const GOBS& g)
{
  _gtype = gobs2gobstype( g );
  _gband = gobs2gobsband( g );
  _gattr = gobs2gobsattr( g );

  return 1;
}

// set from string
// -----------
int t_gobs::gobs(const string& s)
{
  GOBS g = str2gobs(s);  
  _gtype = gobs2gobstype( g );
  _gband = gobs2gobsband( g );
  _gattr = gobs2gobsattr( g );

  return 1;
}

// get GOBS enum
// -----------
GOBS t_gobs::gobs()const
{
  string s = gobstype2str(_gtype) + 
             gobsband2str(_gband) + 
             gobsattr2str(_gattr);
  return str2gobs( s );
}

// get 2char gobs
// ----------
GOBS t_gobs::gobs2CH(GSYS gs)const
{
  GOBS g = X;
  
  if(_gattr == ATTR_NULL) {    // already 2char signal
    g = this->gobs(); 
  }else{
    if(_gtype == TYPE_C){
      if(_gattr != ATTR_C && (gs == GPS || gs == GLO)) g = tba2gobs(TYPE_P, _gband, ATTR_NULL);
      else  g = tba2gobs(_gtype, _gband, ATTR_NULL);
    }    
  }

  return g;
}

// get 2char gobs
// ----------  
GOBS t_gobs::gobs3CH()const
{
  GOBS g = X;
    
  if(_gattr != ATTR_NULL && _gattr != ATTR){   // already 3char signal
    g = this->gobs();
  }else{
         if(_gtype == TYPE_C){ g = tba2gobs(_gtype, _gband, ATTR_C); }
    else if(_gtype == TYPE_P){ g = tba2gobs(TYPE_C, _gband, ATTR_W); }    
  }

  return g;
}
  
// get true if code observation
// -----------
bool t_gobs::is_code()const
{
  return( _gtype == TYPE_C || _gtype == TYPE_P );
}


// get true if phase observation
// -----------
bool t_gobs::is_phase()const
{
  return( _gtype == TYPE_L );
}
/*
// set attributes priority for 2char Pcode to 3char conversion
// ---------------
void t_gobs::_set_attr_priority()
{
  _attr_priority.push_back(ATTR_A);
  _attr_priority.push_back(ATTR_B);
  _attr_priority.push_back(ATTR_C);
  _attr_priority.push_back(ATTR_D);
  _attr_priority.push_back(ATTR_I);
  _attr_priority.push_back(ATTR_L);
  _attr_priority.push_back(ATTR_M);
  _attr_priority.push_back(ATTR_N);
  _attr_priority.push_back(ATTR_P);
  _attr_priority.push_back(ATTR_Q);
  _attr_priority.push_back(ATTR_S);
  _attr_priority.push_back(ATTR_W);
  _attr_priority.push_back(ATTR_X);
  _attr_priority.push_back(ATTR_Y);
  _attr_priority.push_back(ATTR_Z);
}
*/ 
} // namespace
