
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

#include "../gio/gnote.h"
#include "../gutils/gcommon.h"
#include "../gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// constructor
// ---------
t_gnote::t_gnote(t_note n, string f, string s)
{
  gtrace("t_gnote::construct");
  _stat = n;
  _func = f;
  _text = s;
}


// destructor
// ---------
t_gnote::~t_gnote()
{}


// overloading << operator
// -----------------------------
ostream& operator<<(ostream& os, const t_gnote& n)
{
  os << n.str();
  return os;
}


// overloading == operator
// -----------------------------
bool t_gnote::operator==(const t_gnote& n) const
{
//  boost::mutex::scoped_lock lock(_mutex_triple);
  return ( n.status() == _stat &&
	   n.func()   == _func &&
	   n.text()   == _text   
	 );
}


// overloading < operator
// -----------------------------
bool t_gnote::operator<(const t_gnote& n) const
{
//  boost::mutex::scoped_lock lock(_mutex_triple);
  return ( n.status() < _stat &&
	   n.func()   < _func &&
	   n.text()   < _text
	 );
}


// get string
// ---------
string t_gnote::_str()const
{ 
  string note;
  switch( _stat ){
    case GMESSAGE : note = "Message - "; break;
    case GWARNING : note = "Warning - "; break;
    case GERROR   : note = "Error - "  ; break;
  }

  return note;
}

   
// constructor
// ----------
t_gallnote::t_gallnote( )
{
  gtrace("t_gallnote::constructor");
}


// destructor
// ----------
t_gallnote::~t_gallnote()
{
  gtrace("t_gallnote::destructor");
  this->clear();
}


// clean
// --------------------
void t_gallnote::clear()
{
  gtrace("t_gallnote::clear");

  _gmutex.lock();
  _gnotes.clear();
  _gmutex.unlock();
}


// add note (messages/warning/errors)
// ---------------------
void t_gallnote::mesg(t_note note, string func, string text)
{
  gtrace("gallnote::mesg(note, func, text)");
 
  _gmutex.lock();

  t_gnote gnote(note,func,text);

  // eliminate repeating messages
  bool exist = false;
  for( auto it = _gnotes.begin(); it != _gnotes.end(); ++it ){
    if( *it == gnote ){ exist = true; }
  }
  if( !exist ) _gnotes.push_back(gnote);

  _gmutex.unlock(); return;
}


// get note (messages/warning/errors)
// --------------------
vector<t_gnote> t_gallnote::mesg()
{
  gtrace("gallnote::mesg()");

  return _gnotes;
}


} // namespace