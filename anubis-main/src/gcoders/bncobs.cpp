
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <memory>

#include "../gcoders/bncobs.h"
#include "../gutils/gconst.h"
#include "../gutils/gnss.h"
#include "../gutils/gsys.h"
#include "../gall/gallobj.h"
 
using namespace std;

namespace gnut {

// constructor
t_bncobs::t_bncobs( t_gsetbase* s, string version, int sz )
 : t_gcoder( s, version, sz )
{
  _begepoch = false;
//  _endepoch = false;

  _tt = FIRST_TIME;
}

// ---------
// synchronized real-time observation in BNC structure (BNC:bnchelp.html)
// ---------
int t_bncobs::decode_head(char* buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();

  // no header expected, but fill the buffer
  t_gcoder::_add2buffer(buff, sz);
  _mutex.unlock(); return -1;
}


// ---------
// synchronized real-time observation in BNC structure (BNC:bnchelp.html)
// ---------
int t_bncobs::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();
   
//   if( _log ) _log->comment(2,"bncobs","Start decode_data");
   
  // complete gcoder buffer
  if( t_gcoder::_add2buffer(buff, sz) == 0 ){ _mutex.unlock(); return 0; }

//  cout << " BUFFER : \n" << _buffer << "\n size = " << sz  << " END OF BUFFER \n\n"; cout.flush();

  string line;
  int tmpsize = 0;     // individual reading counter
  int consume = 0;     // total read counter
   
  while( ( tmpsize = t_gcoder::_getline( line, 0 ) ) >= 0 ){
//cout << line << endl;
    if ( line.find(crlf) == string::npos ) return 0; // end of buffer
     
    if( line.compare(crlf) == 0 ){  // empty line
      consume += t_gcoder::_consume(tmpsize);
      continue;
    }

    istringstream ss( line );
     
    string c, site;
    int gw, satnum;
    double sow;
    char satsys;
    
    if (line.compare(0,1,">") == 0){ // read epoch
      ss.clear();
      ss >> c >> gw >> sow;
//cout << line << endl;       
//cout << "gw: " << gw << " sow: " << sow << endl;
      if( !ss.fail() ){
        t_gtime t(t_gtime::GPS);
        t.from_gws(gw, sow);
		_tt = t;
		//cout << _tt.str_ymdhms() << " " << ss.str() << endl;
		consume += t_gcoder::_consume(tmpsize);
		continue; 
      }      
    }    

	if (!_validepo(_tt)) {
		consume += t_gcoder::_consume(tmpsize);
		continue;
	}
    
    ss.clear();
    ss >> site >> satsys >> satnum;     
    if (!ss.fail() && (satsys == 'G' || 
             satsys == 'R' || 
             satsys == 'E' ||
             satsys == 'C' ||
             satsys == 'S' ||
             satsys == 'J') ){
      string prn = t_gsys::eval_sat( int2str(satnum), t_gsys::char2gsys(satsys) );

      shared_ptr<t_gobsgnss> obs = make_shared<t_gobsgnss>(site,prn,_tt);

      map<GOBS,double> gobs;
      map<GOBS,int> glli;              
      
      // filtering REC
      if( _rec.size() == 0 || _rec.find( site ) == _rec.end() ){
        consume += t_gcoder::_consume(tmpsize);
        continue;
      }
       
      // checking if receiver object exists
      _check_recobj(site);
      
      // filtering GNSS systems
      if( _sys.size() != 0 && _sys.find( t_gsys::char2str(satsys) ) == _sys.end() ){
        consume += t_gcoder::_consume(tmpsize);
        continue;
      }       
       
      if (satsys == 'R'){
        int slotnum;
        ss.clear();
        ss >> slotnum;
        obs->channel(slotnum);
      }       

      string obstyp;
      double observ;
      int lli;
      ss.clear();       
      
      while( ss >> obstyp >> observ && !ss.fail() ){
		if (t_gsys::char2gsys(satsys) == BDS && obstyp[1]=='1')  obstyp[1] = '2'; // change B1 -> B2 !!!
        GOBS typ = str2gobs(obstyp);
        gobs[typ] = observ;
        
        if (obstyp.compare(0,1,"L") == 0 ) {
          ss >> lli;
          glli[typ] = lli;
        }     
      }

//      map<GOBS,int>::const_iterator itI = glli.begin();
      //while( itI != glli.end() ){ 
      //  if( itI->second != -1 ) obs->addlli(itI->first, itI->second ); itI++;
      //}       

      map<GOBS,double>::const_iterator itD = gobs.begin();
      while( itD != gobs.end() ){ 
        if( itD->second != 0.0 ){                
          obs->addobs(itD->first, itD->second );
        }
        itD++;
      }
       
#ifdef DEBUG       
      cout << fixed << setprecision(3)
           << " DECODED: " << site 
           << " "          << prn
           << " C1C: " << setw(16) << obs->getobs(C1C)
           << " L1C: " << setw(16) << obs->getobs(L1C)
           << " S1C: " << setw(16) << obs->getobs(S1C)        
           << obs->epoch().str(" %Y-%m-%d %H:%M:%S[%T] %W %w %v %J %s") << "\n";
//    int ooo; cin >> ooo;
#endif
       
      // fill data
      map<string,t_gdata*>::iterator it = _data.begin();
      while( it != _data.end() ){
        if (it->second->id_type() == t_gdata::ALLOBS && _validepo(_tt)){
          ((t_gallobs*)it->second)->addobs( obs );
      
		  if (_log) {
			  ostringstream ltmp;
			  double diff = t_gtime::current_time(t_gtime::GPS) - obs->epoch();
			  ltmp << "add site observations " << site << " " << obs->sat() << " "
				   << " " << obs->epoch().str("%Y-%m-%d %H:%M:%S[%T]") << "Diff " << setprecision(2) << diff;
			  _log->comment(3, "bncobs", ltmp.str());
		  }
		}  
        it++;   
      }       
    }
    cnt++;
              
    consume += t_gcoder::_consume(tmpsize);     

  } // while
  _mutex.unlock(); return consume; 

} //end method


bool t_bncobs::_validepo(const t_gtime& t)
{
  if( t > _end || t < _beg ) return false;
  
  if(_filter_epoch(t)) return true;
  else return false;
}

// Checking if a receiver object exists
// ----------------------------------
void t_bncobs::_check_recobj(string& site)
{
  map<string,t_gdata*>::iterator it = _data.begin();
  for( it = _data.begin(); it != _data.end(); ++it ) {
    if( it->second->id_type()  == t_gdata::ALLOBJ ) {
      t_gallobj* all_obj = (t_gallobj*) it->second;
      shared_ptr<t_gobj> rec = all_obj->obj(site);
      if(rec == nullptr) {
        rec = make_shared<t_grec>();
        rec->id(site);                   
        rec->name(site);
        all_obj->add(rec);
      }
    }
  }
  
}
   
} // namespace
