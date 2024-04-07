
#ifndef GNOTE_H
#define GNOTE_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
 
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements motes (messages,warning,errors)
  Version: $ Rev: $

  2017-08-04 /JD: created
  2018-09-14 /JD: updated

-*/

#include <map>
#include <set>
#include <string>
#include <memory>
#include <stdio.h>
#include <fstream>

#include "../gio/glog.h"
#include "../gutils/gmutex.h"

#define BUF_SIZE 1024

using namespace std;

namespace gnut {

enum t_note{ GERROR, GWARNING, GMESSAGE };

class t_gnote {
   
 public:
   t_gnote(t_note n, string f, string s);
   virtual ~t_gnote();
   
   string str()const{ return _str()+_text; }
   string note()const{ return _str(); }
   string text()const{ return _text; }
   string func()const{ return _func; }
   t_note status()const{ return _stat; }

   bool            operator<(const t_gnote& n) const;
   bool            operator==(const t_gnote& n) const;
   friend ostream& operator<<(ostream& os, const t_gnote& n);

 protected:
   virtual string _str()const;
   
   string  _func;          // note function
   string  _text;          // note text
   t_note  _stat;          // note status

 private:
     
};

   
// container for gallnotes, should not be derived from gdata as others   
// ----------
class t_gallnote {

 public:
           t_gallnote();
  virtual ~t_gallnote();

  void mesg(t_note note, string func, string text);     // set/get notes (messages/warning/errors)
  vector<t_gnote> mesg();
   
  void clear();
   
 protected:

  mutable t_gmutex _gmutex;
  vector<t_gnote> _gnotes;                     // cummulate notes message/warning/error

};

} // namespace

#endif
