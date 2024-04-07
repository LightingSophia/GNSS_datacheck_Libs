
#ifndef GLOG_H
#define GLOG_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
 
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements glog class derived from giof
  Version: $ Rev: $

  2011-10-10 /JD: created

-*/

#include <string>
#include <vector>

#include "../gio/giof.h"

#define CACHE_LINES 300

using namespace std;

namespace gnut {

class t_glog : public t_giof {

 public:
   t_glog( string mask = "" );
   virtual ~t_glog();
   
   void comment(int l, const string& str);
   void comment(int l, const string& ide, const string& str);

   void time_stamp(bool b);
   bool time_stamp() const;

   void cache_size(int i);
   int  cache_size() const;
	
   void verb(int i);
   int  verb() const;

   void clear();

 protected:

   bool            _time;             // time stamp
   int             _verb;             // verbosity
   int             _size;             // cache size
   vector<string>  _cache;            // cache for messages
   t_gmutex        _log_gmutex;       // special mutex for comments
   
#ifdef BMUTEX
   boost::mutex    _log_mutex; 
#endif
   
 private:
     
};

} // namespace

#endif
