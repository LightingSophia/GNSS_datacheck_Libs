
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

#include <stdio.h>
#include <math.h>

#include "../gdata/gobj.h"
#include "../gutils/gtypeconv.h"
#include "../gutils/gsysconv.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gobj::t_gobj()
 :  _id(""),
    _name(""),
    _overwrite(false)
{
//  gtrace("t_gobj::constructor");
  id_type(OBJ);
// NOT NECESSARY 
//  t_gtriple eccnull(0.0, 0.0, 0.0);
//  _mapecc[FIRST_TIME] = eccnull;
  
   _pcvnull = 0;
}

/*
// copy constructor
// ----------
t_gobj::t_gobj(const t_gobj& obj)
 :  _id(obj.id()),
    _name(obj.name())
{
   gtrace("t_gobj::constructor");
  _type = obj.id_type();

  vector<t_gtime> vTIM = obj.ecc_id();
  vector<t_gtime>::iterator itTIM = vTIM.begin();
  while( itTIM != vTIM.end() ){
    this->ecc( *itTIM, obj.ecc( *itTIM ));
    itTIM++;
  }

   vTIM = obj.pcv_id();
  itTIM = vTIM.begin();
  while( itTIM != vTIM.end() ){
    this->pcv( *itTIM, obj.pcv( *itTIM ));
    itTIM++;
  }

   vTIM = obj.ant_id();
  itTIM = vTIM.begin();
  while( itTIM != vTIM.end() ){
    this->ant( *itTIM, obj.ant( *itTIM ));
    itTIM++;
  }
}
*/


// destructor
// ----------
t_gobj::~t_gobj()
{
//  gtrace("t_gobj::destructor");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

//_mapantpcv.clear();
  _mappcv.clear();
//_mapeccxyz.clear();
  _mapecc.clear();   
  _mapant.clear();
  _mapcrd.clear();   
//_params.clear();
  _gmutex.unlock();
}


// set id
// ----------
void t_gobj::id( string str )
{
   gtrace("t_gobj::it(str): "+str);   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _id = str;
  _gmutex.unlock();
}


// get name
// ----------
string t_gobj::id() const
{
   gtrace("t_gobj::id");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _id;
  _gmutex.unlock(); return tmp;
}


// set overwrite
// ----------
void t_gobj::overwrite(bool overwrite)
{
   gtrace("t_gobj::overwrite(bool)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _overwrite = overwrite; 
  _gmutex.unlock();
}


// get overwrite
// ----------
bool t_gobj::overwrite()
{
   gtrace("t_gobj::overwrite");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  _gmutex.unlock();   return _overwrite;
}   


// set name
// ----------
void t_gobj::name( string str )
{
   gtrace("t_gobj::name(str): "+str);
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
 _gmutex.lock();

 _name = str;   
   
 _gmutex.unlock(); return; 
}  


// get name
// ----------
string t_gobj::name() const
{
   gtrace("t_gobj::name");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _name;
  _gmutex.unlock(); return tmp;
}


// set domes
// ----------
void t_gobj::domes( string str )
{
   gtrace("t_gobj::domes(str)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  _domes = str;   
   
  _gmutex.unlock(); return; 
}  


// get DOMES
// ----------
string t_gobj::domes() const
{
   gtrace("t_gobj::domes");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _domes;
  _gmutex.unlock(); return tmp;
}


// set description
// ----------
void t_gobj::desc( string str )
{
   gtrace("t_gobj::desc(str)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  _desc = str;   
   
  _gmutex.unlock(); return; 
}   


// get description
// ----------
string t_gobj::desc() const
{
   gtrace("t_gobj::desc");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _desc;
  _gmutex.unlock(); return tmp;
}


// set ecc offsets w.r.t. center of mass/reference point
// ----------
void t_gobj::ecc(const t_gtriple& ecc, const t_ecc_type& typ, const t_gtime& beg, const t_gtime& end )
{
   gtrace("t_gobj::ecc(...)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _ecc(ecc, typ, beg, end);
  _gmutex.unlock(); return; 
}
   
   
// set ecc offsets w.r.t. center of mass/reference point
// ----------
void t_gobj::_ecc(const t_gtriple& ecc, const t_ecc_type& typ, const t_gtime& beg, const t_gtime& end )
{
  t_gtriple zero(0.0, 0.0, 0.0);
  t_gtriple tmp = ecc;
  t_mapecc::iterator it =  _mapecc.find(beg);

  if( typ == ECC_XYZ ){
  	t_gtriple crd = this->_crd(beg);
    if( crd.zero() ){
      string lg = "Warning: " + _id + " no valid coordinates, cannot store eccentricity !";
      return;
    }
   	t_gtriple ell; xyz2ell(crd, ell, false);
    t_gtriple xyz; neu2xyz(ell, tmp, xyz); 
    tmp = xyz;
  }
  
  if( end < beg ){
     string lg = "Warning: " + _id + " not valid end time (end<beg) for eccentricity neu (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";
     if( _log ) _log->comment(0,"gobj", lg );
     else               cerr << "gobj: " << lg << endl;
     return;
  }
   
  // begin record
  if( it == _mapecc.end() ){     // not exists
      _mapecc[beg] = tmp;
     
  }else{                         // record exists
    if( it->first == LAST_TIME ||
        it->second == zero ){

       _mapecc[beg] = tmp;
    }else{
       string lg = "Warning: " + _id + " valid object record cannot be changed for eccentricity xyz";
       if( _log && !_overwrite ) _log->comment(1,"gobj", lg );
       return;
    }
  }

  // control end of record (with new beg search)
  it = _mapecc.find(beg); it++;
   
  // beg was last in map (add final empty record)
  if( it == _mapecc.end() ){
    _mapecc[end] = zero;

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
        if( it->second == zero ) _mapecc.erase(it); // remove obsolete empty record
        _mapecc[end] = zero;
      
      }else{                                        // too close to next record
        string lg = "Warning: object " + _id + " ECC end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
        else               cerr << "gobj: " << lg << endl;
      }
    }else if (end != it->first){
      string lg = "Warning: object " + _id + " ECC end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
      if( _log ) _log->comment(2,"gobj",lg );
    }
  }

  // remove duplicated empty records
  auto itNEW = _mapecc.begin();
  auto itOLD = itNEW;
  while( itOLD != _mapecc.end() ){    
    if( ++itNEW != _mapecc.end() ){
      if( ( itNEW->second == zero && itOLD->second == zero ) ){
        _mapecc.erase( itNEW++ );
      }
    }
    itOLD = itNEW;
  }
  return;
}


// get ecc offsets (>=t) w.r.t. center of mass/reference point (interface only)
// ----------
t_gtriple t_gobj::ecc(const t_gtime& t, const t_ecc_type& typ) const
{
   gtrace("t_gobj::ecc");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple ecc(0.0, 0.0, 0.0);
  ecc = this->_ecc(t, typ);

  _gmutex.unlock(); return ecc;
}
   

// get ecc offsets (>=t) w.r.t. center of mass/reference point
// ---------
t_gtriple t_gobj::_ecc(const t_gtime& t, const t_ecc_type& typ) const
{
	t_gtriple ecc(0.0, 0.0, 0.0);
	t_gtriple neu(0.0, 0.0, 0.0);
	t_mapecc::const_iterator it = _mapecc.upper_bound(t);

	if( it != _mapecc.begin()){
    it--; ecc = neu = it->second;
  }

  if( typ == ECC_XYZ && !neu.zero() ){
	  t_gtriple crd = this->_crd(t);
		if( !crd.zero()){
     	t_gtriple ell;
      xyz2ell(crd, ell, false);
			neu2xyz(ell, neu, ecc);
		}
  }
	return ecc;
}


// return validity for eccneu at epoch t
// ----------
void t_gobj::ecc_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const
{
  gtrace("t_gobj::ecc_validity");   
  
  t_mapecc::const_iterator it = _mapecc.upper_bound(t);
  if( it != _mapecc.begin() && it != _mapecc.end()){
    end = it->first;
    it--;
    beg = it->first;
  }
}


// set object position
// ----------
void t_gobj::crd(const t_gtriple& crd, const t_gtriple& std, const t_gtime& beg, const t_gtime& end )
{
  gtrace("t_gobj::crd(t_gtriple, epo)");

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  _crd(crd, std, beg, end);
   
  _gmutex.unlock(); return; 
}


// set object position
// ----------
void t_gobj::_crd(const t_gtriple& crd, const t_gtriple& std, const t_gtime& beg, const t_gtime& end )
{
  pair<t_gtriple, t_gtriple> val;
  val = make_pair(crd, std);
  
  pair<t_gtriple, t_gtriple> zero;
  t_gtriple nullCRD(0.0, 0.0, 0.0);
  t_gtriple nullSTD(0.0, 0.0, 0.0);
  zero = make_pair(nullCRD, nullSTD);
  
  t_mapcrd::iterator it =  _mapcrd.find(beg);
   if( end < beg ){   
      string lg = "Warning: " + _id + " not valid end time (end<beg) for coordinates (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";
      if( _log ) _log->comment(0,"gobj",lg );
      else               cerr << "gobj: " << lg << endl;
      return;
  }
  // begin record
  if( it == _mapcrd.end() ){     // not exists
      _mapcrd[beg] = val;
  }else{                         // record exists
    if( it->first == LAST_TIME ||
        it->second == zero ){
       _mapcrd[beg] = val;
    }else{
       string lg = "Warning: " + _id + " valid object record cannot be changed (coordinates)";
       if( _log && !_overwrite ) _log->comment(1,"gobj",lg );
       return;
    }
  }
  // control end of record (with new beg search)
  it = _mapcrd.find(beg); it++;
   
  // beg was last in map (add final empty record)
  if( it == _mapcrd.end() ){
    _mapcrd[end] = zero;

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
        if( it->second == zero ) _mapcrd.erase(it); // remove obsolete empty record
        _mapcrd[end] = zero;
        
      }else{                                        // too close to next record
        string lg = "Warning: object " + _id + " 'crd' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
        else               cerr << "gobj: " << lg << endl;
      }
    }else if(end != it->first){
      string lg = "Warning: object " + _id + " 'crd' end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
      if( _log ) _log->comment(2,"gobj",lg );
    }
  }

  // remove duplicated empty records
  t_mapcrd::iterator itNEW = _mapcrd.begin();
  t_mapcrd::iterator itOLD = itNEW;
  while( itOLD != _mapcrd.end() ){    
    if( ++itNEW != _mapcrd.end() ){
      if( ( itNEW->second == zero && itOLD->second == zero ) )
        //          ( itNEW->first == LAST_TIME ) )
      {
        _mapcrd.erase( itNEW++ );
      }
    }
    itOLD = itNEW;
  }

  return;
}


// get fixed position (using eccentricity)
// ----------
t_gtriple t_gobj::crd_fix(const t_gtime& t) const
{
	gtrace("t_gobj::crd_fix");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	// applying eccentricity
	t_gtriple crd(0.0, 0.0, 0.0);	
	double dt = 0.0;
	t_mapcrd mapcrd;
	for (auto iter = _mapcrd.begin(); iter != _mapcrd.end();iter++ ) {
		t_gtriple std = iter->second.second;
		if (0.1 < std.norm() && std.norm() < 999 || iter->second.first.norm() == 0.0) 			
			continue;		
		if (dt == 0.0 || fabs(iter->first - t) < dt) {
			dt = fabs(iter->first - t);
			crd = iter->second.first;
			mapcrd[iter->first]=iter->second;
		}
	}

	t_gtriple marker(0.0, 0.0, 0.0);
	marker = this->_crd(t);

	// applying eccentricity
	t_gtriple arp(0.0, 0.0, 0.0);
	arp = crd + this->_ecc(t, ECC_XYZ);

	_gmutex.unlock(); return arp;
}


// get object position (using eccentricity)
// ----------
t_gtriple t_gobj::crd_arp(const t_gtime& t) const
{
//   gtrace("t_gobj::crd_arp");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple marker(0.0, 0.0, 0.0);
  marker = this->_crd(t);

  // applying eccentricity
  t_gtriple arp(0.0, 0.0, 0.0);
  arp = marker + this->_ecc(t,ECC_XYZ);
   
  _gmutex.unlock(); return arp;
}


// interface for _crd(const t_gtime& t) const
// ----------
t_gtriple t_gobj::crd(const t_gtime& t) const
{
   gtrace("t_gobj::crd");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple marker(0.0, 0.0, 0.0);
  marker = this->_crd(t);

  _gmutex.unlock(); return marker;
}


// interface for _std(const t_gtime& t) const
// ----------
t_gtriple t_gobj::std(const t_gtime& t) const
{
   gtrace("t_gobj::std");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple marker(0.0, 0.0, 0.0);
  marker = this->_std(t);

  _gmutex.unlock(); return marker;
}


// get object position (without eccentricity correction)
// ----------
t_gtriple t_gobj::_crd(const t_gtime& t) const
{  
#ifdef DEBUG
   for( t_mapcrd::const_iterator itE = _mapcrd.begin(); itE != _mapcrd.end(); ++itE  )
     cerr << "OBJ-CRD: searching EPOCH: " << t.str_ymdhms() 
          << "   found: " << itE->first.str_ymdhms() << "\n";
#endif

  t_gtriple crd(0.0, 0.0, 0.0);
  t_mapcrd::const_iterator it = _mapcrd.upper_bound(t);
  if( it == _mapcrd.begin() ){ return crd; }   // not found
  it--;
  t_gtriple tmp = it->second.first;

  return tmp;
}


// get object position std
// ----------
t_gtriple t_gobj::_std(const t_gtime& t) const
{  
#ifdef DEBUG
   for( t_mapcrd::const_iterator itE = _mapcrd.begin(); itE != _mapcrd.end(); ++itE  )
     cerr << "OBJ-CRD: searching EPOCH: " << t.str_ymdhms() 
          << "   found: " << itE->first.str_ymdhms() << "\n";
#endif

  t_gtriple std(0.0, 0.0, 0.0);
  t_mapcrd::const_iterator it = _mapcrd.upper_bound(t);
  if( it == _mapcrd.begin() ){ return std; }   // not found
  it--;
  t_gtriple tmp = it->second.second;

  return tmp;
}   


// return validity for crd at epoch t
// ----------
void t_gobj::crd_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const
{
  gtrace("t_gobj::crd_validity");   

#ifdef DEBUG  
  for( t_mapcrd::const_iterator itE = _mapcrd.begin(); itE != _mapcrd.end(); ++itE  )
    cerr << "OBJ-CRD: searching EPOCH: " << t.str_ymdhms()
        << "         testing EPOCH: "  << itE->first.str_ymdhms()
        << fixed << setprecision(3) << setw(16) << (itE->second).first[0] << endl;
#endif

  t_mapcrd::const_iterator it = _mapcrd.upper_bound(t);
  if( it != _mapcrd.begin() ){
    end = it->first;
    it--;
    beg = it->first;
  }
}


// return epochs
// ----------
set<t_gtime> t_gobj::crd_epochs() const
{
 _gmutex.lock();
 set<t_gtime> tmp;
 for( auto it = _mapcrd.begin(); it != _mapcrd.end(); ++it  ){
   tmp.insert(it->first);
 }
 _gmutex.unlock();
 return tmp;
}

// set pcv element
// ----------
void t_gobj::pcv(shared_ptr<t_gpcv> pcv, const t_gtime& beg, const t_gtime& end)
{
   gtrace("t_gobj::pcv(...)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  this->_pcv(pcv, beg, end);
  _gmutex.unlock();
}


// set pcv element (protected)
// ----------
void t_gobj::_pcv(shared_ptr<t_gpcv> pcv, const t_gtime& beg, const t_gtime& end)
{
   gtrace("t_gobj::_pcv");   

  // ------------------------------------------------------------------
  // ------------------------------------------------------------------
  // _mappcv[beg] = pcv;  // _mappcv[t] = pcv->pointer(); !!!! ???????? 
  // ------------------------------------------------------------------
  // ------------------------------------------------------------------
  
  /* t_mapantpcv::iterator itant = _mapantpcv.find(ant); */
  /* if (itant == _mapantpcv.end()){ */
  /*    _mapantpcv[ant][beg] = pcv; */
  /* }else{ */
  /*   t_mappcv mappcv = itant->second; */
     
    t_mappcv::iterator it = _mappcv.find(beg);

    if( end < beg ){
      string lg = "Warning: " + _id + " not valid end time (end<beg) for PCV (beg:" + beg.str_ymdhms() + " -> " + end.str_ymdhms() + ")";
      if( _log ) _log->comment(0,"gobj",lg );
      else               cerr << "gobj: " << lg << endl;
      return;
    }
   
    // begin record
    if( it == _mappcv.end() ){     // not exists
      _mappcv[beg] = pcv;        // ?? POINTER ??
    
    // record exists
    }else if( it->first == LAST_TIME || !it->second ){
      _mappcv[beg] = pcv;        // ?? POINTER ??

    }else{
      string lg = "Warning: " + _id + " valid object record cannot be changed (PCV)";
      if( _log && !_overwrite ) _log->comment(1,"gobj",lg );
      return;
    }

    // control end of record (with new beg search)
    it = _mappcv.find(beg); it++;
     
    // beg was last in map (add final empty record)
    if( it == _mappcv.end() ){
      _mappcv[end] = _pcvnull;
    }else{                                            // process end according to next value
      if( end < it->first ){                          // only if end is smaller then existing
        if( fabs(it->first - end) > 3600 ){           // significantly smaller!
          if( it->second == 0 ) _mappcv.erase(it); // remove obsolete empty record
          _mappcv[end] = _pcvnull;
          
        }else{                                        // too close to next record
          string lg = "Warning: object " + _id + " 'pcv' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
          if( _log ) _log->comment(2,"gobj",lg );
        }
      }else if(end != it->first){
        string lg = "Warning: object " + _id + " 'pcv' end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
      }
    }

    // remove duplicated empty records
    t_mappcv::iterator itNEW = _mappcv.begin();
    t_mappcv::iterator itOLD = itNEW;
    while( itOLD != _mappcv.end() ){    
      if( ++itNEW != _mappcv.end() ){
        if( ( !itNEW->second && !itOLD->second ) )
          //            ( itNEW->first == LAST_TIME ) )
        {
          _mappcv.erase( itNEW++ );
        }
      }
      itOLD = itNEW;
    }     
//  }
}


// get pcv element (>=t)
// ----------
shared_ptr<t_gpcv> t_gobj::pcv(const t_gtime& t) const
{  
   gtrace("t_gobj::pcv(t)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  shared_ptr<t_gpcv> tmp = _pcv(t);
  _gmutex.unlock(); return tmp;
}


// get pcv element (>=t) (protected)
// ----------
shared_ptr<t_gpcv> t_gobj::_pcv(const t_gtime& t) const
{  
   gtrace("t_gobj::_pcv(t)");   

  shared_ptr<t_gpcv> pcv;   
  string ant;      

  t_mapant::const_iterator it = _mapant.upper_bound(t);
  if( it == _mapant.begin() ) {
    ant = "";  // antenna not found
    string lg = "Warning: unknown PCO (no antenna found in the object " + _id + " ) " + t.str_ymdhms();
    if( _log ) _log->comment(1,"gobj",lg );
    return _pcvnull;
  }else{
    ant = (--it)->second;     
  }

  t_mappcv::const_iterator it2 = _mappcv.upper_bound(t);
  if( it2 == _mappcv.begin() ) {     
     ostringstream ostr;
     ostr << "Warning: unknown PCO ( antenna " << left << setw(20) << ant << " not found in ATX ) " << t.str_ymdhms();
     if( _log ) _log->comment(1,"gobj",ostr.str() );
     return _pcvnull;
  }else{ 
    pcv = (--it2)->second;
  }

  if(pcv->pcvkey().compare(ant) != 0){
     if(pcv->pcvkey().compare(0, 16, ant, 0, 16) == 0){
       string lg = "Warning: PCO Used without considering randome " + pcv->pcvkey() + " " + ant;
       if (_log) _log->comment(1, "gobj", lg);
     }
     else{
       ostringstream ostr;
       ostr << "Warning: unknown PCO ( changed antenna " << left << setw(20) << ant << " not found in ATX ) " << t.str_ymdhms();
       if (_log) _log->comment(1, "gobj", ostr.str());
#ifdef DEBUG     
       cout << "gobj RESET: ant [" << ant << "] [" << pcv->pcvkey() << " " << t.str_ymdhms(" epo:") << "]\n";
#endif     
       return _pcvnull;
    }
  }   
  return pcv;
}


// set antenna name
// ----------
void t_gobj::ant(string ant, const t_gtime& beg, const t_gtime& end)
{
   gtrace("t_gobj::and(...)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  _ant(ant, beg, end); 
   
  _gmutex.unlock(); return; 
}
   

// set antenna name
// ----------
void t_gobj::_ant(string ant, const t_gtime& beg, const t_gtime& end)
{  
  t_mapant::iterator it =  _mapant.find(beg);

  if( end < beg ){
    string lg = "Warning: " + _id + " not valid end time (end<beg) for antenna (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";
    if( _log ) _log->comment(0,"gobj",lg );
    else               cerr << "gobj: " << lg << endl;
    return;
  }

  // begin record
  if( it == _mapant.end() ){     // not exists
      _mapant[beg] = ant;
     
  }else{                         // last value
    if( it->first == LAST_TIME ||
        it->second.empty() ){

       _mapant[beg] = ant;
    }else{                      // record exists
       string lg = "Warning: " + _id + " valid object record cannot be changed (antenna)";
       if( _log && !_overwrite ) _log->comment(1,"gobj",lg );
       return;
    }
  }

  // control end of record (with new beg search)
  it = _mapant.find(beg); it++;

  // beg was last in map (add final empty record)
  if( it == _mapant.end() ){
    _mapant[end] = "";

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
        if( it->second.empty() ) { 
          _mapant.erase(it); // remove obsolete empty record
          _mapant[end] = "";
        }else {
          _mapant[end] = it->second;
        }
        
      }else{                                        // too close to next record
        string lg = "Warning: object " + _id + " 'obj' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
        else               cerr << "gobj: " << lg << endl;
      }
    }else if(end != it->first){
      string lg = "Warning: object " + _id + " 'ant' " + ant + " end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
      if( _log ) _log->comment(2,"gobj",lg );
    }
  }

  // remove duplicated empty records
  t_mapant::iterator itNEW = _mapant.begin();
  t_mapant::iterator itOLD = itNEW;
  while( itOLD != _mapant.end() ){    
    if( ++itNEW != _mapant.end() ){
      if( ( itNEW->second.empty() && itOLD->second.empty() ) )
        //          ( itNEW->first == LAST_TIME ) )
      {
        _mapant.erase( itNEW++ );
      }
    }
    itOLD = itNEW;
  }

  return;
}


// get antenna name (>=t)
// ----------
string t_gobj::ant(const t_gtime& t) const
{  
   gtrace("t_gobj::ant()");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = this->_ant(t);
  _gmutex.unlock(); return tmp;
}


// get antenna name (>=t) (protected)
// ----------
string t_gobj::_ant(const t_gtime& t) const
{  
   gtrace("t_gobj::_ant(t)");   
   
  t_mapant::const_iterator it = _mapant.upper_bound(t);
  if( it == _mapant.begin() ) return "";  // not found
  it--;
  return it->second;
}

     
// get time tags
// ----------
/*vector<t_gtime> t_gobj::ecc_id() const
{
  vector<t_gtime> tmp;
  t_mapecc_xyz::const_iterator itMAP = _mapeccxyz.begin();
  while( itMAP != _mapecc.end() ){
    tmp.push_back( itMAP->first );
    itMAP++;
  }
  return tmp;
}
*/


// return validity for antenna at epoch t
void t_gobj::ant_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const
{
  gtrace("t_gobj::ant_validity");   
  
  t_mapant::const_iterator it = _mapant.upper_bound(t);
  if( it != _mapant.begin() ){
    end = it->first;
    it--;
    beg = it->first;
  }
}


// get time tags
// ----------
vector<t_gtime> t_gobj::pcv_id() const
{      
   gtrace("t_gobj::pcv_id");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<t_gtime> tmp;
//  t_mapantpcv::const_iterator itant = _mapantpcv.find(ant);
   
//  if (itant == _mapantpcv.end()){
//     return tmp;
//  }else{
  t_mappcv::const_iterator itMAP = _mappcv.begin();
  while( itMAP != _mappcv.end() ){
    tmp.push_back( itMAP->first );
    itMAP++;
  }
  _gmutex.unlock(); return tmp;     
//  }   
}


// get time tags
// ----------
vector<t_gtime> t_gobj::ant_id() const
{
   gtrace("t_gobj::ant_id");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  vector<t_gtime> tmp = this->_ant_id();
  _gmutex.unlock(); return tmp;   
}


// get time tags (protected)
// ----------
vector<t_gtime> t_gobj::_ant_id() const
{   
   gtrace("t_gobj::_ant_id");   
   
  vector<t_gtime> tmp;
  t_mapant::const_iterator itMAP = _mapant.begin();
  while( itMAP != _mapant.end() ){
    tmp.push_back( itMAP->first );
    itMAP++;
  }
  return tmp;
}


// set mappcv for all t_gobj
// -----------------------
void t_gobj::sync_pcv(t_gallpcv* pcvs)
{
  gtrace("t_gobj::sync_pcv");   
   
  if( pcvs == 0 ) return;

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<t_gtime> vant = this->_ant_id();

  if( vant.size() == 0 ){  // ant in obj is not set
    if( _log ) _log->comment(1,"gobj","Warning: " + _id + " ant name not set: " );
    else                                    cerr << "gobj - Warning: " + _id + " ant name not set:" << endl;

  }else{
    for( vector<t_gtime>::iterator it = vant.begin(); it != vant.end(); it++ ){

      // set pcv for all antennas
      string ant = this->_ant(*it);
      t_gtime epo = *it;  
   
      shared_ptr<t_gpcv> gpcv = pcvs->gpcv(ant, "*", epo); 
      if( gpcv != _pcvnull ) {
        this->_pcv(gpcv, epo);
      }else{
        ostringstream ostr;
        ostr << "Warning: unknown PCO ( antenna " << left << setw(20) << ant << " not found in ATX ) " << epo.str_ymdhms();
        if( _log && !ant.empty() ) _log->comment(1,"gobj",ostr.str() );
        else if( _log &&  ant.empty() ) _log->comment(2,"gobj",ostr.str() );    
      }     
   
#ifdef DEBUG
       cout << " gobj.cpp - sync for ant " << ant << "  pcv found [" << gpcv << "] " << it->str_ymdhms() << endl;
#endif
    }
  }

  _gmutex.unlock(); return;
}


// check consistency
// ----------
void t_gobj::compare(shared_ptr<t_gobj> gobj, const t_gtime& tt, string source)
{
  gtrace("t_gobj::compare");
  
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif

  _gmutex.lock();
  
  string txt = " "+_id+"/"+source+" ";

  string old, alt;
  t_gtriple trip_old, trip_alt;
  t_gtriple std_old, std_alt;  
  t_gtime beg, end;

  // Name
  old = trim(_name); alt = trim(gobj->name());

  if( old != alt && !alt.empty()){
    if( old.empty() ){
      _name = alt;
      if( _log ) _log->comment(1,"gobj","Object completed," + txt + tt.str_ymdhms("NAME: ") + " -> " + alt);
    }else if( _overwrite ){
      _name = alt;
      if( _log ) _log->comment(1,"gobj","Object updated," + txt + tt.str_ymdhms("NAME: ") + " -> " + alt);
    }else{
      if( _log ) _log->comment(1,"gobj","Object not match," + txt + tt.str_ymdhms("NAME: ") + " -> " + old + " vs " + alt);
    }
  }
   
  // Domes
  old = trim(_domes); alt = trim(gobj->domes());
  
  if( old != alt && !alt.empty()){
    if( old.empty() ){
      _domes = alt;
      if( _log ) _log->comment(1,"gobj","Object completed," + txt + tt.str_ymdhms("DOME: ") + " -> " + alt);
    }else if( _overwrite ){
      _domes = alt;
      if( _log ) _log->comment(1,"gobj","Object updated," + txt + tt.str_ymdhms("DOME: ") + " -> " + alt);
    }else{
      if( _log ) _log->comment(1,"gobj","Object not match," + txt + tt.str_ymdhms("DOME: ") + " -> " + old + " vs " + alt);
    }
  }   
   
  // Antenna
  old = trim(_ant(tt)); alt = trim(gobj->ant(tt));
  gobj->ant_validity(tt, beg, end);
  
  if( old != alt && !alt.empty()){
    if( old.empty() ){
      _ant(alt, beg, end);
      this->ant_validity(tt, beg, end);
      if( _log ) _log->comment(1,"gobj","Object completed," + txt + beg.str_ymdhms("ANT: ") + end.str_ymdhms(" -> ") + " -> " + alt);
    }else if( _overwrite ) {
      _ant(alt, beg, end);
      if( _log ) _log->comment(1,"gobj","Object updated," + txt + beg.str_ymdhms("ANT: ") + end.str_ymdhms(" -> ") + " -> " + alt);
    }else{
      if( _log ) _log->comment(1,"gobj","Object not match," + txt + tt.str_ymdhms("ANT: ") + " -> " + old + " vs " + alt);
    }
  }

  // Coordinates
  trip_old = _crd(tt); trip_alt = gobj->crd(tt);
  std_old  = _std(tt); std_alt  = gobj->std(tt);
  gobj->crd_validity(tt, beg, end);

  t_gtriple diff = trip_old - trip_alt;
  double dist = diff.norm();
   
  if( (dist > 10 || std_alt < std_old) && !trip_alt.zero() ) {   // fixed max diff 50 m
    if( trip_old.zero() ){
      _crd(trip_alt, std_alt, beg, end);
      this->crd_validity(tt, beg, end);
      if( _log ) _log->comment(1,"gobj","Object completed," + txt + beg.str_ymdhms("CRD: ") + end.str_ymdhms(" -> ")
                                         + dbl2str(trip_alt[0]) + dbl2str(trip_alt[1]) + dbl2str(trip_alt[2]));
    }else if( _overwrite ){
      _crd(trip_alt, std_alt, beg, end);
      if( _log ) _log->comment(1,"gobj","Object updated," + txt + beg.str_ymdhms("CRD: ") + end.str_ymdhms(" -> ")
                                         + dbl2str(trip_alt[0]) + dbl2str(trip_alt[1]) + dbl2str(trip_alt[2]));
    }else{
      if( _log ) _log->comment(1,"gobj","Object not match," + txt + tt.str_ymdhms("CRD: "));
    }
  }

  // Eccentricity   
  trip_old = _ecc(tt,ECC_NEU); trip_alt = gobj->ecc(tt,ECC_NEU);
  gobj->ecc_validity(tt, beg, end);
  
  if( trip_old != trip_alt && !trip_alt.zero()){
    if( trip_old.zero() ){
      _ecc(trip_alt, ECC_NEU, beg, end);
      this->ecc_validity(tt, beg, end);
      if( _log ) _log->comment(1,"gobj", "Object completed," + txt + beg.str_ymdhms("ECC: ") + end.str_ymdhms(" -> ")
                                         + dbl2str(trip_alt[0],4) + dbl2str(trip_alt[1],4) + dbl2str(trip_alt[2],4));
    }else if( _overwrite ){
      _ecc(trip_alt, ECC_NEU, beg, end);
      if( _log ) _log->comment(1,"gobj","Object updated," + txt + beg.str_ymdhms("ECC: ") + end.str_ymdhms(" -> ")
                                         + dbl2str(trip_alt[0],4) + dbl2str(trip_alt[1],4) + dbl2str(trip_alt[2],4));
    }else{
      if( _log ) _log->comment(1,"gobj","Object not match," + txt + tt.str_ymdhms("ECC: "));
    }
  }
/*
  // ECC_XYZ eccentricity
  trip_old = _ecc(tt,ECC_XYZ); trip_alt = gobj->ecc(tt,ECC_XYZ);
  gobj->eccxyz_validity(tt, beg, end);
    
  if( trip_old != trip_alt && !trip_alt.zero()){
    if ( trip_old.zero() ) {
      _ecc(trip_alt, ECC_XYZ, beg, end);
      this->eccxyz_validity(tt, beg, end);
      if( _log ) _log->comment(1,"gobj","Object completed," + txt + beg.str_ymdhms("ECC: ") + end.str_ymdhms(" -> ")
                                         + dbl2str(trip_alt[0],4) + dbl2str(trip_alt[1],4) + dbl2str(trip_alt[2],4));
    }else if( _overwrite ){     
      _ecc(trip_alt, ECC_XYZ, beg, end);
      if( _log ) _log->comment(1,"gobj","Object updated," + txt + beg.str_ymdhms("ECC: ") + end.str_ymdhms(" -> ")
                                         + dbl2str(trip_alt[0],4) + dbl2str(trip_alt[1],4) + dbl2str(trip_alt[2],4));
    }else{
      if( _log ) _log->comment(1,"gobj","Object not match," + txt + tt.str_ymdhms("ECC: "));
    }
  }   
*/
  _gmutex.unlock(); return;
}


// operator for sorting
// ----------
bool t_gobj::operator<(const t_gobj& t) const 
{ 
 return _name < t.name();
}
   

// operator for sorting
// ----------
bool t_gobj::operator==(const t_gobj& t) const
{
 return  _name == t.name();
}
   
} // namespace
