
#ifndef GTRN_H
#define GTRN_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements transmitter class
  Version: $ Rev: $

  2013-06-20 /PV: created
  2018-08-05 /JD: updated

-*/

#include <stdio.h>
#include <string>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "../gdata/gobj.h"
#include "../gdata/gdata.h"
#include "../gdata/grxnhdr.h"

using namespace std;

namespace gnut {

class t_gtrn : public t_gobj {

 public:
   t_gtrn();
//   t_gtrn(const t_gtrn& obj);
   virtual ~t_gtrn();

//   typedef map<t_gtime, int>        t_mapchk;  // for changing channel canal of Glonass

   typedef pair<string,t_rxnhdr>            t_header_pair;
   typedef vector<t_header_pair>            t_header;
      
   virtual void     header(const t_rxnhdr& hdr, string path);
   virtual t_rxnhdr header(string path) const;
   virtual t_header headers() const;

   virtual bool isrec() {return true;}
   virtual bool istrn() {return false;}
   
//   void chk(int chk, const t_gtime& beg, const t_gtime& end);
//   int  chk(const t_gtime& t) const;       // set/get receiver
  
   virtual void channel(int chk);
   virtual int  channel() const;
   
   
 protected:
   t_rxnhdr        _header(const t_gtime& epo);

   t_header        _headers;               // map of rinex header information
//   t_mapchk        _mapchk;
  int _channel;   // temporarily only one value. Must be enhance via _mapchk 
 private:
};

} // namespace

#endif
